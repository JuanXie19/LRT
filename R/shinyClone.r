#' Explanatory analysis of clonotype
#'
#' @param Combined a \code{CombinedDataSet} object
#' @import Seurat
#' @import dplyr
#' @import DT
#' @import ggplot2
#' @import shiny
#'
#' @return shiny app
#'
#' @export


shinyClone <- function(CombinedDataSet =NULL){

# extract needed information
	META <- CombinedDataSet@seurat@meta.data
	NAMES <- colnames(META)
	LEN <- function(x){
		if(class(x)=='character' & length(unique(x))<=11){return(1)}
		else{return(0)}
		}
	tt <- apply(META,2,LEN)
	GROUPVAR <- NAMES[which(tt==1)]  ## find the variables that could serve grouping varaible

	cloneCall <- 'CTaa'
	colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

	### define argument used in repertoire analysis
	
	ui <-navbarPage('shinyClone',
	
		 tabPanel('Clonotype',fluidRow(
			column(6,
				selectInput(inputId = 'shownTable',
							label = 'show distribution of clonal type across',
							choices = GROUPVAR,
							selected = 'clusters'),
				hr(),
				DT::dataTableOutput('clonetable')),
			column(6,
				radioButtons(inputId = 'visualize',
								label = 'Visualize',
								choices = c('UMAP'='umap','t-SNE'='tsne','PCA'='pca'),
								selected = 'umap',inline=T),
				hr(),
								
				fluidRow(
				column(6,
				#h4('Dimensional reduction plot'),
				plotOutput('ReduceDimPlot'),
				textOutput(outputId = 'desc'),
				downloadButton('downloadPlot',label = 'Download plot')),
				column(6,
				#h4('Density contour plot'),
				plotOutput('DensityPlot'),
				textOutput(outputId = 'densitydesc'),
				uiOutput("mybutton")
				)
				)))),
		 tabPanel('Repertoire',fluidRow(
			column(3,
				wellPanel(
					selectInput(inputId = 'analysisType',
								label = 'Choose analysis',
								choices = c('Occupied repertoire','Chord diagram','Diversity','Proportion','Overlap'),
								selected = 'Occupied repertoire'),
					hr(),
					selectInput(inputId = 'groupVariable',
									label = 'Group clonal type by',
									choices = GROUPVAR,
									selected = 'clusters'),
					hr(),
					bsCollapse(id = 'analysisArgument',open = 'Parameters',
						bsCollapsePanel('Parameters',
							uiOutput('repertorieDetails')
							)))),								

			column(9,
				plotOutput('repertoirePlot'),
				textOutput('text1'),
				downloadButton('downloadRepertoireplot',label = 'Download plot')
				))
		 
		 ))



	server <-function(input,output,session){
		
		
		CTaa.df <- reactive({
				tmp <- META %>% select(CTaa, group = input$shownTable)
				tmp1 <- data.frame(rbind(table(tmp$CTaa,tmp$group)))
				return(tmp1)
		})
		cell.num <- reactive(rowSums(CTaa.df()))
		group.num <- reactive(apply(CTaa.df(),1, function(x) sum(x>0)))
			
		

		output$clonetable = DT::renderDataTable(data.frame(CTaa.df(),cell.num=cell.num(),group.num=group.num()),server=FALSE,extensions = 'Buttons',options = list(
										dom = 'frtipB',buttons = list('copy','print',list(extend = 'collection',buttons = c('csv','excel'),text='Download table'))))

  
		s <- reactive({
			req(input$clonetable_rows_selected)
			input$clonetable_rows_selected
			})
	
	
########## clonal type analysis
	output$desc <-renderText({	
		CLONE <- rownames(CTaa.df())[s()]		
		paste('This clonotype consists of ',cell.num()[s()],'cells', 'spanning',group.num()[s()],'group/groups')
	})
	
  plotUMAP <- function(){
  
	
	PLOT <-reactive({
	
		REDUCED <- CombinedDataSet@seurat@reductions[[input$visualize]]@cell.embeddings
		CTaa.df <- META%>%select(CTaa,group = input$shownTable)
		DATA <- data.frame(CTaa.df,Dim1 = REDUCED[,1],Dim2=REDUCED[,2])
		INDEX <- which(DATA$CTaa == rownames(CTaa.df())[s()])
		G.df.highlight <- DATA[INDEX,]
		
		switch(input$visualize,
		'umap' = ggplot(DATA,aes(Dim1,Dim2,color=group))+geom_point(alpha=0.3,shape=1,size=0.5)+
				geom_point(data=G.df.highlight,aes(Dim1,Dim2),size=1.2,shape=17)+theme_bw()+ xlab("UMAP_1") + ylab("UMAP_2")+
				xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1),max(DATA$Dim1)), ylim = c(min(DATA$Dim2),max(DATA$Dim2)), expand = FALSE)+
			theme(legend.position="bottom"),
		'tsne' = ggplot(DATA,aes(Dim1,Dim2,color=group))+geom_point(alpha=0.3,shape=1,size=0.5)+
				geom_point(data=G.df.highlight,aes(Dim1,Dim2),size=1.2,shape=17)+theme_bw()+ xlab("tSNE_1") + ylab("tSNE_2")+
				xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1),max(DATA$Dim1)), ylim = c(min(DATA$Dim2),max(DATA$Dim2)), expand = FALSE)+
			theme(legend.position="bottom"),
		'pca' = ggplot(DATA,aes(Dim1,Dim2,color=group))+geom_point(alpha=0.3,shape=1,size=0.5)+
				geom_point(data=G.df.highlight,aes(Dim1,Dim2),size=1.2,shape=17)+theme_bw()+ xlab("PC_1") + ylab("PC_2")+
				xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1),max(DATA$Dim1)), ylim = c(min(DATA$Dim2),max(DATA$Dim2)), expand = FALSE)+
			theme(legend.position="bottom")
		)
	})
	return(PLOT())
  }

	
    output$ReduceDimPlot <-renderPlot({
		if(length(s())){
		print(plotUMAP())
	  }
    })
	
	output$downloadPlot <-downloadHandler(
		filename = function(){
		CLONE <- rownames(CTaa.df())[s()]
		paste(CLONE,'_reducedDimension.pdf',sep='')},
		content = function(file){
			 pdf(file)
			 print(plotUMAP())
			 dev.off()	
	})
	
	
	DENSITY <- function(){
	
	KDE2D <-reactive({
		REDUCED <- CombinedDataSet@seurat@reductions[[input$visualize]]@cell.embeddings
		CTaa.df <- META%>%select(CTaa,group = input$shownTable)
		DATA <- data.frame(CTaa.df,Dim1 = REDUCED[,1],Dim2=REDUCED[,2])
		INDEX <- which(DATA$CTaa == rownames(CTaa.df())[s()])
		G.df.highlight <- DATA[INDEX,]
	
		if(nrow(G.df.highlight)>=10){
			switch(input$visualize,
			'umap' = ggplot(G.df.highlight,aes(Dim1,Dim2))+geom_point(size=0.2)+geom_density_2d() + xlab("UMAP_1") + ylab("UMAP_2")+
			theme_bw()+ xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1),max(DATA$Dim1)), ylim = c(min(DATA$Dim2),max(DATA$Dim2)), expand = FALSE),
			'tsne' = ggplot(G.df.highlight,aes(Dim1,Dim2))+geom_point(size=0.2)+geom_density_2d() + xlab("tSNE_1") + ylab("tSNE_2")+
			theme_bw()+ xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5), ylim = c(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5), expand = FALSE),
			'pca' = ggplot(G.df.highlight,aes(Dim1,Dim2))+geom_point(size = 0.2)+geom_density_2d() + xlab("PC_1") + ylab("PC_2")+
			theme_bw()+ xlim(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5)+ylim(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5)+
			coord_cartesian(xlim = c(min(DATA$Dim1)-0.5,max(DATA$Dim1)+0.5), ylim = c(min(DATA$Dim2)-0.5,max(DATA$Dim2)+0.5), expand = FALSE)
			)
		}
	})
	return(KDE2D())
	
	}
	
	
	output$DensityPlot <-renderPlot({
		if(length(s())){
		print(DENSITY())
	  }
    })
	
	output$densitydesc <- renderText({	
		REDUCED <- CombinedDataSet@seurat@reductions[[input$visualize]]@cell.embeddings
		CTaa.df <- META%>%select(CTaa,group = input$shownTable)
		DATA <- data.frame(CTaa.df,Dim1 = REDUCED[,1],Dim2=REDUCED[,2])
		INDEX <- which(DATA$CTaa == rownames(CTaa.df())[s()])
		
		if(length(INDEX)>=8){
			paste('The plot shows the density contour of cells. Dots denote cells, and curves denote contour.')
		}
		else{
			paste('')
			}
		
	})
	
	output$mybutton <- renderUI({
		REDUCED <- CombinedDataSet@seurat@reductions[[input$visualize]]@cell.embeddings
		CTaa.df <- META%>%select(CTaa,group = input$shownTable)
		DATA <- data.frame(CTaa.df,Dim1 = REDUCED[,1],Dim2=REDUCED[,2])
		INDEX <- which(DATA$CTaa == rownames(CTaa.df())[s()])
		
		if(length(INDEX)>=10){
			downloadButton('downloadDensity', 'Download plot')
		}
	})
	
	output$downloadDensity <-downloadHandler(
		filename = function(){
		CLONE <- rownames(CTaa.df())[s()]
		paste(CLONE,'_densityplot.pdf',sep='')},
		content = function(file){
			 pdf(file)
			 print(DENSITY())
			 dev.off()	
		})
	
############   repertoire analysis
 
	
  output$repertorieDetails <- renderUI({
			if(is.null(input$analysisType))
				return()
			
			switch(
				input$analysisType,
				'Overlap' = selectInput(inputId = 'overlapIndex',
										label = 'method',
										choices = c("overlap", "morisita", "jaccard", "raw"),
										selected = 'morisita'),
				
				'Occupied repertoire' = list(radioButtons(inputId = 'showFreq',
												label = 'show option',
												choices = c('show frequency','show proportion'),
												selected = 'show frequency'),
				
											numericInput(inputId = 'clone_Small',
															label = 'Small',
															min =2,max =10,value=5),
											numericInput(inputId = 'clone_Medium',
															label = 'Medium',
															min = 11,max = 50 ,value = 20),
											numericInput(inputId = 'clone_Large',
															label = 'Large',
															min = 51, max = 150, value = 100),
											numericInput(inputId = 'clone_Hyperexpand',
														label = 'Hyperexpand',
														min = 151,max = 600,value = 500)
											),
				'Chord diagram' = radioButtons(inputId = 'showFreq2',
												label = 'show option',
												choices = c('show frequency','show proportion'),
												selected = 'show frequency')
										)			
		})
		
		observeEvent(input$showFreq=='show proportion',{
			updateNumericInput(session,'clone_Small', min =0.001,max =0.009,value=0.001,step=0.001)
			updateNumericInput(session,'clone_Medium',min = 0.01,max = 0.09 ,value = 0.01,step = 0.01)
			updateNumericInput(session,'clone_Large',min = 0.1, max = 0.5, value = 0.1,step = 0.1)
			updateNumericInput(session,'clone_Hyperexpand',min = 0.6, max = 1, value = 1,step = 0.1)
		})
		
		observeEvent(input$showFreq=='show frequency',{
			updateNumericInput(session,'clone_Small', min =2,max =10,value=5,step=1)
			updateNumericInput(session,'clone_Medium',min = 11,max = 50 ,value = 20,step = 1)
			updateNumericInput(session,'clone_Large',min = 51, max = 150, value = 100)
			updateNumericInput(session,'clone_Hyperexpand',min = 151, max = 600, value = 500)
		})
	
	METHOD <- reactive(input$overlapIndex)

	group.by <- reactive(input$groupVariable)
	
	PROP <- reactive({
				switch(input$showFreq,
				'show proportion' = TRUE,
				'show frequency' = FALSE)})
	
	PROP2 <- reactive({
				switch(input$showFreq2,
				'show proportion' = TRUE,
				'show frequency' = FALSE)})
	
	clone_small <- reactive(input$clone_Small)
	clone_medium <- reactive(input$clone_Medium)
	clone_large <- reactive(input$clone_Large)
	clone_hyperexpand <- reactive(input$clone_Hyperexpand)
	
	CIRCLE <- function(){

		temp <- reactive({
			circle_temp <- LRT::getCirclize(Combined,group.by=group.by(),proportion=PROP2())
			grid.cols <- scales::hue_pal()(length(unique(Combined@seurat@active.ident)))
			names(grid.cols) <- levels(Combined@seurat@active.ident)
			circlize::chordDiagram(circle_temp,self.link = 1, grid.col = grid.cols)
		})
		return(temp())
		}
		
	REPERTOIRE <- function(){
		TEMP <- reactive({
			switch(input$analysisType,
				'Occupied repertoire' = LRT::occupiedscRepertoire(Combined,group.by = group.by(),proportion=PROP(),cloneTypes = c(Small=clone_small(),Medium=clone_medium(),Large=clone_large(),Hyperexpanded = clone_hyperexpand()),label = TRUE,facet.by = NULL, na.include = FALSE,exportTable = FALSE) ,
				'Chord diagram' = print(CIRCLE()),
				'Diversity' = LRT::clonalDiversity(Combined,group.by = group.by(), x.axis = NULL, split.by = NULL,exportTable = FALSE, n.boots = 100),
				'Proportion' = LRT::clonalProportion(Combined,group.by = group.by(),split = c(10, 100, 1000, 10000, 50000), 
                      exportTable = FALSE),
				'Overlap' = LRT::clonalOverlap(Combined,group.by = group.by(),method = METHOD(), exportTable = FALSE)
			)
		})
		return(TEMP())
	}
	
	output$text1 <- renderText({
		switch(EXPR = paste(input$analysisType,input$showFreq2),
		'Chord diagram show frequency' = paste('This plot visualizes the interconnection of groups, i.e., the number of clonotypes shared between groups'),
		'Chord diagram show proportion' = paste('This plot visualizes the interconnection of groups, i.e., the proportion of clonotypes shared between groups')
		)
	
	})
	output$repertoirePlot <- renderPlot({print(REPERTOIRE())})
	output$downloadRepertoireplot <- downloadHandler(
		filename = function(){
				paste(input$analysisType,'_plot.pdf',sep='')
			},
			content = function(file){
				pdf(file)
				print(REPERTOIRE())
				dev.off()			
			}
	)

}

## create shiny object
shinyApp(ui = ui, server = server)
}