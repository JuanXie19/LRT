#' cluster inferred clonotype trajectories
#' @param Combined a \code{CombinedDataSet} object
#' @param trajectory a \code{TrajectoryDataSet} object
#'
#' @export
#' @import dplyr
#' @import shiny
#' @import shinyBS
#' @import shinysky
#' @import ggplot2
#' @import parallel
#' @import DirichletMultinomial
#' @import tidyr
#' @import slingshot
#' @import Seurat
#' @import SingleCellExperiment
#'
#' @examples
#' TCR <-read.csv("/PATH/TO/YOUR/scTCR-seqData/",header=T)
#' load('Mice.sub2.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getOverallTrajectory(Combined,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE)
#' shinyClust(Trajectory)




shinyClust <- function(TrajectoryDataSet = NULL){
	library(SingleCellExperiment)
	
	PT <- as.data.frame(TrajectoryDataSet@slingPseudotime) 
	LIN.name <- colnames(PT)[grep('Lin',colnames(PT))]
	
	ui <- fluidPage(
	 # App title ----
	titlePanel("shinyClust"),

	fluidRow(
	 column(2,
		wellPanel(	
			actionButton('btn_go','Run clustering'),
			hr(),
			bsCollapse(id="clusterInfo", open='Clustering options',

                   bsCollapsePanel("Clustering options",
								sliderInput(inputId = 'clustMax',
											label = 'Maximum number of clonotype clusters',
											min = 2, max = 10, value = 7),
								radioButtons(inputId = 'dmmCluster',
										label = 'Way to decide number of clusters',
										choices = c('data driven','user specified'),
										selected = 'data driven'),
								uiOutput('clusteringDetails'),
								radioButtons(inputId = "criteria",
                                               label = "Criteria to check model fit",
                                               choices = c("Laplace","BIC", "AIC"),
                                               selected = 'Laplace')
															
								)                                  
                   )					   			
				)),
		column(10,
			
		tabsetPanel(
			tabPanel('Clustering',
			fluidRow(
				busyIndicator(text = "Running clonotype clustering, which may take a few minutes, please wait..."),
				column(6,
					uiOutput('fitPlot'),
					uiOutput('fitText')

				),
				column(6,
					DT::dataTableOutput('posteriorTable'),
					uiOutput('posteriorText')
				))),
			tabPanel('Clonotype cluster exploration',
			fluidRow(
				selectInput(inputId = 'clusterIndex',
							label = 'Choose a clonotype cluster to display',
							choices = c('clonotype cluster 1','clonotype cluster 2' ),
							selected = 'clonotype cluster 1')
				),
			fluidRow(
				column(6,
				h4('Density contour'),
				uiOutput(outputId = 'plotDen'),
				uiOutput('downloadDen'),
				helpText('This plot shows the density contour for a clonotype cluster,overlaid with overall trajectory.',
					'Blue curves represent density contour.', 'Orange/green curve represent overall trajectory')),

				column(6,
				h4('Top 10 clones'),
				uiOutput(outputId = 'plotBar'),
				uiOutput('downloadBar'),
				helpText('This plot shows the phenotypic distribution of cells from the top 10 expanded clones.', 
				'Colors denote cell type/state.'))
		
		
			)
		),
		tabPanel('Biasedness evaluation',
			fluidRow(
				fluidRow(
					column(6,selectInput(inputId = 'lineageIndex',
							label = 'Choose a lineage to display',
							choices = LIN.name,
							selected = 'Lineage 1')),
					column(6,sliderInput(inputId = 'permNum',,
									label = 'number of permutations',
									min = 20, max = 10000,
									value = 5000))),				
				fluidRow(
					column(6,plotOutput('plotPseudo'),
					helpText('The plot shows the distribution of pseudotimes.')),
					column(6,
					busyIndicator(text = "Performing permutation test, which may take a few minutes, please wait..."),
					DT::dataTableOutput('permuTable'),
					helpText('The table shows the p value for permutation tests.')
					)
					)
			)),
		tabPanel('Clonotype cluster characterization',
			fluidRow(
				column(5,
					plotOutput('diversityPlot'),
					helpText('The plot shows the diversity indices for each clonotype cluster')
				),
				column(7,
					radioButtons(inputId = "geneInd",
                                  label = "Show",
                                               choices = c("V gene","J gene"),
                                               selected = 'V gene',inline=T),
					plotOutput('geneUsage'),
					helpText('The plot shows the gene usage for each clonotype cluster')
				))
		)	
		)))
	

	)
	
	server <- function(input,output,session){

		output$clusteringDetails <- renderUI({
			switch(input$dmmCluster,
				'user specified' = numericInput(inputId = 'clustNumS',
												label = 'Desired number of clonotype clusters',
												min = 2, value =6)
			)
		})
		
		########## clustering
		seurat <- TrajectoryDataSet@seurat
		temp.df <- data.frame(cell = colnames(seurat),clone = seurat$CTaa,cell.cluster = seurat$clusters)
		temp.df2 <- temp.df %>% dplyr::group_by(clone,cell.cluster) %>% dplyr::summarise(n = n())
		temp.df2.wide <- temp.df2 %>% pivot_wider(names_from = cell.cluster,values_from = n,values_fill=0)
		temp.df2.wide2 <- temp.df2.wide[rowSums(temp.df2.wide[,2:ncol(temp.df2.wide)])>4,]
		
		
		kD <- reactive(input$clustMax)
		kS <- reactive(input$clustNumS)
		

################		
#### run clonotype clustering
	
		myDMM <- eventReactive(input$btn_go,{
				mclapply(1:kD(), dmn,count = as.matrix(temp.df2.wide2[,-1],verbose = TRUE))		
			})
		
		

		observeEvent(myDMM(),{
			updateSelectInput(session,'clusterIndex',choices = paste('clonotype cluster',1:length(myDMM())))
		})
		
		
		
		clusterRST <- eventReactive(input$btn_go,{
				if (input$dmmCluster =='data driven')
					x <- switch(input$criteria,
						'Laplace' = {
									lplc <- sapply(myDMM(),laplace)
									best <- myDMM()[[which.min(lplc)]]
									LABELS <- mixture(best,assign = T)
									post <- round(fitted(best,scale = T),3)
									colnames(post) <- paste0('clonotype cluster',1:ncol(post))
									list(LABELS,post)
									},
						'AIC' = {
									aic <- sapply(myDMM(),AIC)
									best <- myDMM()[[which.min(aic)]]
									LABELS <- mixture(best,assign = T)
									post <- round(fitted(best,scale = T),3)
									colnames(post) <- paste0('clonotype cluster',1:ncol(post))
									list(LABELS,post)},
						'BIC' = {
									bic <- sapply(myDMM(),BIC)
									best <- myDMM()[[which.min(bic)]]
									LABELS <- mixture(best,assign = T)
									post <- round(fitted(best,scale = T),3)
									colnames(post) <- paste0('clonotype cluster',1:ncol(post))
									list(LABELS,post)})	
				if (input$dmmCluster == 'user specified')
					x <- list(mixture(myDMM()[[kS()]],assign = T),fitted(myDMM()[[kS()]],scale = T))
				return(x)
						
		})
		
		
### clustering plot
	observeEvent(input$btn_go,{
			output$fitPlot	<- renderUI({
				switch(input$criteria,
					'Laplace' = plotOutput('laplace'),
					'AIC' = plotOutput('aicPlot'),
					'BIC' = plotOutput('bicPlot')
				)
			})		
		})
		
	output$laplace <- renderPlot({
			lplc <- sapply(myDMM(),laplace)
			plot(lplc,type = 'b',xlab = 'Number of clonotype clusters', ylab = 'Laplace criterion value')
		})
		
	output$aicPlot <- renderPlot({
			aic <- sapply(myDMM(),AIC)
			plot(aic,type = 'b',xlab = 'Number of clonotype clusters', ylab = 'AIC')
		})
		
	output$bicPlot <- renderPlot({
			bic <- sapply(myDMM(),BIC)
			plot(bic,type = 'b',xlab = 'Number of clonotype clusters', ylab = 'BIC')
		})
		
	observeEvent(input$btn_go,{
		output$posteriorTable = DT::renderDataTable(data.frame(clusterRST()[[2]]))	
		})		
	
	## description of plot and table
	output$fitPlotText <- renderText({
			paste("The plot shows how model fit changes with the number of clonotype clusters. We'd prefer a model with lower model fit value.")
		})
	
	observeEvent(input$btn_go,{
		output$fitText = renderUI({
				textOutput('fitPlotText')
		   })	
		})
	
	output$postTableText <- renderText({
			paste("The table shows the mean posterior probabilities of cells from a clonotype cluster belong to a cell state.")
		})
	
	observeEvent(input$btn_go,{
		output$posteriorText = renderUI({
				textOutput('postTableText')
		   })	
		})

################# cluster exploration 
	# update the number of clusters based on clustering results
	observeEvent(clusterRST(),{
			updateSelectInput(session,'clusterIndex',choices = paste('clonotype cluster',1:max(clusterRST()[[1]])))
		})
	clusterID <- reactive(as.numeric(strsplit(input$clusterIndex, ' ')[[1]][3]))  # used on cluster exploration
	
	
## density contour
	REDUCED <- seurat@reductions[['umap']]@cell.embeddings
	DATA <- data.frame(CTaa = seurat$CTaa,group = seurat$clusters,Dim1 = REDUCED[,1],Dim2=REDUCED[,2])
	 
	 
	dims = seq_len(2)
	linInd <- NULL
	if(is.null(linInd)){
		linInd <- seq_along(TrajectoryDataSet@slingLineages)
	}

	DATA.list <- list()
	for(ii in seq_along(TrajectoryDataSet@slingCurves)[linInd]){
		c <- TrajectoryDataSet@slingCurves[[ii]]
		DATA.list[[ii]] <- c$s[c$ord,dims]
 
  #lines(c$s[c$ord, dims], lwd = lwd, col = col[ii])
	}
	names(DATA.list) <- paste0('curve',seq_along(TrajectoryDataSet@slingCurves)[linInd])

	DATA.df <- as.data.frame(do.call(rbind,DATA.list))
	DATA.df$curves <- rep(seq_along(TrajectoryDataSet@slingCurves)[linInd],c(sapply(DATA.list,nrow)))


	DEN <- function(CLONES){
		INDEX <- which(DATA$CTaa %in% CLONES)
		G.df.highlight <- DATA[INDEX,]
		ggplot(G.df.highlight,aes(Dim1,Dim2))+geom_density_2d(linewidth=0.9) +
			geom_path(data = DATA.df,aes(UMAP_1,UMAP_2,color=factor(curves)),linewidth = 1)+ 
			xlab("UMAP_1") + ylab("UMAP_2")+
			xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+
			ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+theme_classic()	
	}
	
	observeEvent(input$btn_go,{
			output$plotDen <- renderUI({
				plotOutput('densityPlot')
			})
			
		})
		
		
		
	output$densityPlot <- renderPlot({
		IND <- which(clusterRST()[[1]]==clusterID())
		CLONES <- temp.df2.wide2[IND,]$clone
		print(DEN(CLONES))
		})
		
		## download density plots
		observeEvent(input$btn_go,{
			output$downloadDen <-renderUI({
				downloadButton('downloadDensity',label = 'Download plot')
			})
		})
		

		
		output$downloadDensity <- downloadHandler(
			filename = function(){
				paste(input$clusterIndex,'_density.pdf',sep='')},
				content = function(file){
					pdf(file)
					IND <- which(clusterRST()[[1]]==clusterID())
					CLONES <- temp.df2.wide2[IND,]$clone
					print(DEN(CLONES))
					dev.off()
				}			
		)
		

## top10 clones


		
	
	observeEvent(input$btn_go,{
			output$plotBar <- renderUI({
				plotOutput('barPlot')
			})
			
		})
		
				
	output$barPlot <- renderPlot({
		cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
		temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
		temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]	
		a  <- temp.df.update %>% dplyr::group_by(clone.cluster,clone) %>% dplyr::summarise(n = n()) %>% dplyr::arrange(desc(n),.by_group = T) %>% top_n(10)
		DIST <- list()
		for (i in 1:nrow(a)){
			clone <- a[i,]$clone
			DIST[[i]] <- temp.df2.wide2[which(temp.df2.wide2$clone==clone),2:ncol(temp.df2.wide2)]
		}

		DIST <- do.call(rbind.data.frame,DIST)
		a <- cbind(a,DIST)
		a.sub <- a[which(a$clone.cluster == clusterID()),]
		a.sub.long <- melt(a.sub[,-3],id.vars = c('clone.cluster','clone'))
		a.sub.long$clone = factor(a.sub.long$clone, levels = a.sub$clone)


		ggplot(a.sub.long, aes(fill=variable, y=value, x=clone)) + 
			geom_bar(position="stack", stat="identity")+
			xlab('Top 10 clones')+ylab('Clone size')+theme_classic()+ 
			theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
		
				
		})

## download bar plots
		observeEvent(input$btn_go,{
			output$downloadBar <-renderUI({
				downloadButton('downloadBarplot',label = 'Download plot')
			})
		})
		

		
		output$downloadBarplot <- downloadHandler(
			filename = function(){
				paste(input$clusterIndex,'_bar.pdf',sep='')},
				content = function(file){
					pdf(file)
					cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
					temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
					temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]	
					a  <- temp.df.update %>% dplyr::group_by(clone.cluster,clone) %>% dplyr::summarise(n = n()) %>% dplyr::arrange(desc(n),.by_group = T) %>% top_n(10)
					DIST <- list()
					for (i in 1:nrow(a)){
						clone <- a[i,]$clone
						DIST[[i]] <- temp.df2.wide2[which(temp.df2.wide2$clone==clone),2:ncol(temp.df2.wide2)]
					}

					DIST <- do.call(rbind.data.frame,DIST)
					a <- cbind(a,DIST)
					a.sub <- a[which(a$clone.cluster == clusterID()),]
					a.sub.long <- reshape2::melt(a.sub[,-3],id.vars = c('clone.cluster','clone'))
					a.sub.long$clone = factor(a.sub.long$clone, levels = a.sub$clone)


				ggplot(a.sub.long, aes(fill=variable, y=value, x=clone)) + 
					geom_bar(position="stack", stat="identity")+
					xlab('Top 10 clones')+ylab('Clone size')+theme_classic()+ 
					theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())	
					
					dev.off()
				}			
		)
		
		
#######################
###### baisedness evaluation
    ## distribution of pseudotim
		

		output$plotPseudo <- renderPlot({
			cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
			temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
			temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]	
			
			
			PT$cloneCluster <- factor(temp.df.update$clone.cluster)
			
			ggplot(PT, aes(x = Lineage1, fill = stat(x))) +
				geom_histogram() +facet_grid(rows = vars(cloneCluster))+
				scale_fill_viridis_c(name = "Pseudotime")+
				theme_classic()
				
		})
		
	### bias evaluation
		P <- reactive(input$permNum)
	    
		linID <- reactive(as.numeric(gsub("\\D", "", input$lineageIndex))) 
		lineage <- reactive(input$lineageIndex)
		
		LineageNum <- length(LIN.name)  # check how many lineages it has
		output$permuTable <- DT::renderDataTable(DT::datatable({
		
			cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
			temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
			temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]	
			PT$CLUSTERS <- seurat$clusters
			PT$cloneCluster <- factor(temp.df.update$clone.cluster)
		
			PT.sub <- PT%>%select(lineage(), CLUSTERS,cloneCluster)
			PT.sub <- PT.sub[!is.na(PT.sub[,1]),]   # filter NA rows
			obs.median <- aggregate(PT.sub[,1],list(PT.sub$cloneCluster),FUN = median)
			obs.IQR <- aggregate(PT.sub[,1],list(PT.sub$cloneCluster),FUN = IQR)
			
			
		    set.seed(123456)
			rst_median <- vector('list',P())
		    rst_IQR <- vector('list',P())
			
			for (i in 1:P()){
				x <- split(sample(PT.sub[,1]),rep(1:length(unique(PT.sub$cloneCluster)),c(unclass(table(PT.sub$cloneCluster)))))
				rst_median[[i]] <- unlist(lapply(x,function(s) median(s,na.rm = T))) 
				rst_IQR[[i]] <- unlist(lapply(x,function(s) IQR(s,na.rm=T)))  
			}

			RST_median <- do.call(rbind.data.frame,rst_median)
			RST_IQR <- do.call(rbind.data.frame,rst_IQR)
			
			p.IQR <- vector()
			for (i in 1:nrow(obs.IQR)){
				p.IQR[i] <- sum(RST_IQR[,i] <obs.IQR[i,2])/P()
			
			}
			
			p.median <- vector()
			for (i in 1:nrow(obs.median)){
				MIN <- min(RST_median[,i])
				MAX <- max(RST_median[,i])
				dist1 <- abs(MIN-obs.median[i,2])
				dist2 <- abs(MAX-obs.median[i,2])
				p.median[i] <- ifelse(dist1 < dist2, sum(RST_median[,i] < obs.median[i,2])/P(),sum(RST_median[,i] > obs.median[i,2])/P())
			
			}
			
			### lineage bias
			## only if there are more than 1 lineage
			if(length(LIN.name)>1){
				perm.data <- data.frame(clusters = seurat$clusters,clone.cluster = temp.df.update$clone.cluster)
				
				IND <- which(!is.na(PT[,1]) &!is.na(PT[,2]))   # the cells shared by both lineages
				
				## this is just a compromise, ideally should calculate the ratio of cells in two lineages
				LARGE1 <- head(PT[order(-PT[,1]),], 100)
				LARGE2 <- head(PT[order(-PT[,2]),], 100)
				
				A <- names(which.max(unclass(table(LARGE1$CLUSTERS))))
				B <- names(which.max(unclass(table(LARGE2$CLUSTERS))))
				
				cal_ratio <- function(x){
					log2(length(which(x==A))/length(which(x==B)))
				}
				
				set.seed(123456)
				results <- vector('list',P())
				for (i in 1:P()){
					x <- split(sample(perm.data$clusters),rep(1:length(unique(PT$cloneCluster)),c(unclass(table(temp.df.update$clone.cluster)))))
					results[[i]]<- unlist(lapply(x,cal_ratio))
				}

				RST <- do.call(rbind.data.frame,results)

				RATIO <- unclass(table(PT$cloneCluster,PT$CLUSTERS))
				obs.ratio <- log2(RATIO[,A]/RATIO[,B]) 
				p.Lineage <- vector()
				for (i in 1:length(obs.ratio)){
				MIN <- min(RST[,i])
				MAX <- max(RST[,i])
				dist1 <- abs(MIN-obs.ratio[i])
				dist2 <- abs(MAX-obs.ratio[i])
				p.Lineage[i] <- ifelse(dist1 < dist2, sum(RST[,i] < obs.ratio[i])/P(),sum(RST[,i] > obs.ratio[i])/P())			
				}
				bias.RST <- data.frame(p.Lineage,p.IQR,p.median)
				rownames(bias.RST) <- paste0('clonotype cluster',1:nrow(obs.IQR))
				colnames(bias.RST) <- c('Lineage','Dynamic','Early/Late')		
			}else{
				bias.RST <- data.frame(p.IQR,p.median)
				rownames(bias.RST) <- paste0('clonotype cluster',1:nrow(obs.IQR))
				colnames(bias.RST) <- c('Dynamic','Early/Late')					
			}
			bias.RST				
		}))
		
################ repertoire


######### Diversity
	

	output$diversityPlot <- renderPlot({
		group.by = 'clusters'
		split.by = NULL
		x.axis = 'clone.cluster'
		n.boots = 300
		cloneCall <- "CTaa"
		colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
	
		cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
		temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
		temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]	
		seurat$clone.cluster <- temp.df.update$clone.cluster
		df <- expression2List_seurat(seurat,group.by = 'clusters')
		cloneCall <- theCall(cloneCall)
		df <- checkBlanks(df)
 
 
		mat <- NULL
		mat_a <- NULL
		sample <- c()
		if (!is.null(group.by) || !is.null(x.axis)) {
		df <- bind_rows(df, .id = "element.names")
		df$group.element <- paste0(df[,group.by], ".", df[,x.axis])
		#group.element.uniq <- unique(df$group.element)
		df <- split(df, f = df[,"group.element"])
		}
		min <- short.check(df)
		for (i in seq_along(df)) {
			data <- as.data.frame(table(df[[i]][,cloneCall]))
			mat_a <- NULL
		sample <- c()
   
			for (j in seq(seq_len(n.boots))) {
				x <- sample_n(data, min)
				sample <- diversityCall(x)
				mat_a <- rbind(mat_a, sample)
			}
		mat_a[is.na(mat_a)] <- 0
		mat_a<- colMeans(mat_a)
		mat_a<-as.data.frame(t(mat_a))
		mat <- rbind(mat, mat_a)
		}
	colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE", "Inv.Pielou")
	mat[,"Inv.Pielou"] <- 1 - mat[,"Inv.Pielou"]
	if (!is.null(group.by)) {
	mat[,group.by] <- stringr::str_split(names(df), "[.]", simplify = TRUE)[,1]
	} else {
	group.by <- "Group"
	mat[,group.by] <- names(df)
	}
	if (!is.null(x.axis)) {
		mat[,x.axis] <- stringr::str_split(names(df), "[.]", simplify = TRUE)[,2]
	} else {
		x.axis <- "x.axis"
		mat[,x.axis] <- 1
	}
	rownames(mat) <- names(df)
 
	melt <- suppressMessages(melt(mat, id.vars = c(group.by, x.axis)))
	
	
		ggplot(melt, aes(x=clone.cluster, y=as.numeric(value)))+
		geom_boxplot(outlier.alpha = 0) +
			geom_jitter(aes(color = clusters), size = 3) + 
			labs(color="Group") +
			ylab("Index Score") +
			facet_wrap(~variable, scales = "free", ncol = 2) +
			theme_classic() + 
			theme(axis.title.x = element_blank())
		
	})


########### gene usage
	#gene = reactive(strsplit(input$geneInd,' '))[[1]][1]
	gene = 'V'
   chain = "TRB" 
   plot = "heatmap" 
   plot = 'bar'
   y.axis = "clone.cluster" 
   order = "gene"
   scale = TRUE 
   group.by = 'clone.cluster'
   split.by = NULL
   exportTable = FALSE   

   
	

	
	output$geneUsage <- renderPlot({
		
		cloneCluster.rst <- data.frame(clone=temp.df2.wide2$clone,clone.cluster = clusterRST()[[1]])
		temp.df.update <- left_join(temp.df,cloneCluster.rst,by='clone')
		temp.df.update <- temp.df.update[temp.df.update$clone%in% temp.df2.wide2$clone,]
		seurat$clone.cluster <- temp.df.update$clone.cluster

		df <- expression2List_seurat(seurat,group.by = 'clone.cluster')
		if(!is.null(group.by)) {
			df <- groupList(df, group.by)
		}
		for(i in seq_along(df)) {
			df[[i]] <- off.the.chain(df[[i]], chain, "CTaa")
		}
		df <- bound.input.return(df)
		if (y.axis %!in% colnames(df) | is.null(y.axis)) {
			if (y.axis %!in% c("V", "D", "J", "C")) {
			y.axis <- "element.names"
		} else {
			df <- select.gene(df, chain, y.axis)
			colnames(df)[ncol(df)] <- y.axis
			}
		}
		df <- select.gene(df, chain, gene)
		df <- subset(df, !is.na(df[,ncol(df)])) #remove NA values
		df <- subset(df, df[,ncol(df)] != "NA") #remove values that are character "NA"
		df <- subset(df, df[,ncol(df)] != "") #remove rows with non genes
   #df <- table(df[,ncol(df)], df[,y.axis])
   
		if (!is.null(y.axis) && y.axis != "element.names") {
			df <- df %>%
			dplyr::group_by(df[,ncol(df)], df[,y.axis], element.names) %>%
			dplyr::count() 
		} else {
			df <- df %>%
			dplyr::group_by(df[,ncol(df)], element.names) %>%
			dplyr::count() 
		}
		df <- df %>% dplyr::group_by(element.names) %>% dplyr::mutate(sum = sum(n))
		col.lab <- "Total n"
		if (scale == TRUE) {
			df[,"n"] <- df[,"n"]/df[,"sum"]
			col.lab <- "Scaled Values"
		} 
		colnames(df)[1:2] <- c("Var1", "Var2")
		df <- df %>% dplyr::group_by(Var1, Var2) %>% dplyr::mutate(varcount = sum(n),sd = sd(n, na.rm = TRUE),mean = mean(n))
		if (order == "variance") {
			varOrder <- order(df$varcount, decreasing = TRUE)
			df$Var1 <- factor(df$Var1, levels = unique(df$Var1[varOrder]))
		}
   
	
	
		df2 <- unique(df[,c("Var1", "Var2", "sd", "mean")])
		ggplot(df2, aes(x=Var1, y = mean)) + 
		geom_bar(stat = "identity") + 
		geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                     position=position_dodge(.9)) + 
		theme_classic() + 
		theme(axis.title.x = element_blank(), #remove titles
             axis.title.y = element_blank(), 
             axis.ticks.x = element_blank(), #removes ticks
             axis.text.x = element_text(angle = 90, 
                                        vjust = 0.5, hjust=1, size=rel(0.8))) + 
       facet_grid(Var2~.)

	})

		
	}
	
	shinyApp(ui = ui, server = server)

}