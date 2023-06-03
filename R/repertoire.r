#' repertorie analysis functions used for shinyClone 
#' modified based on the functions from Borcherding N's scRepertoire R package

#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of the 
#' overlap coefficient, morisita, or jaccard index for clonotypes. 
#' The overlap coefficient is calculated using the intersection of clonotypes 
#' divided by the length of the smallest component. 
#'
#' @examples
#' 
#' Combined <- getCombinedDataSet(TCR,Mice.sub2)
#' clonalOverlap(combined,method = "overlap",exportTable = FALSE)
#'
#' @param Combined The \code{CombinedDataSet} object obtained from \code{getCombinedDataSet} function.
#' @param method The method to calculate the overlap, either the "overlap" 
#' coefficient, "morisita", "jaccard" indices, or "raw" for the base numbers.
#' @param exportTable Returns the data frame used for forming the graph
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list

clonalOverlap <- function(Combined,group.by, method = c("overlap", "morisita", "jaccard", "raw"),  
                                exportTable = FALSE){
	cloneCall <- "CTaa"
	colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

    df <- list.input.return(Combined,group.by)
    #cloneCall <- theCall(cloneCall)
    #df <- checkBlanks(df, cloneCall)
    df <- df[order(names(df))]
    values <- stringr::str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]
    num_samples <- length(df[])
    names_samples <- names(df)
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples
    length <- seq_len(num_samples)
    
    if (method == "overlap") {
        coef_matrix <- overlapIndex(df, length, coef_matrix)
    } else if (method == "morisita") {
        coef_matrix <- morisitaIndex(df, length, coef_matrix)
    } else if (method == "jaccard") {
        coef_matrix <- jaccardIndex(df, length, coef_matrix)
    } else if (method == "raw") {
        coef_matrix <- rawIndex(df, length, coef_matrix)
    }
    coef_matrix$names <- rownames(coef_matrix)
    if (exportTable == TRUE) { return(coef_matrix) }
    coef_matrix <- suppressMessages(reshape2::melt(coef_matrix))[,-1]
	
    col <- colorblind_vector(7)
    coef_matrix$variable <- factor(coef_matrix$variable, levels = values)
    coef_matrix$names <- factor(coef_matrix$names, levels = values)
    plot <- ggplot(coef_matrix, aes(x=names, y=variable, fill=value)) +
            geom_tile() + labs(fill = method) +
            geom_text(aes(label = round(value, digits = 3))) +
            scale_fill_gradientn(colors = rev(colorblind_vector(5)), na.value = "white") +
            theme_classic() + theme(axis.title = element_blank())
    return(plot) }
	
#' Visualize the number of single cells with clonotype frequencies by cluster
#'
#' View the count of clonotypes frequency group in seurat or SCE object 
#' meta data after combineExpression(). The visualization will take the 
#' new meta data variable "cloneType" and plot the number of cells with
#' each designation using a secondary variable, like cluster. Credit to 
#' the idea goes to Drs. Carmona and Andreatta and their work with
#' \href{https://github.com/carmonalab/ProjecTILs}{ProjectTIL}.
#'
#' @examples
#' Combined <- getCombinedDataSet(TCR, Mice.sub2)
#' occupiedscRepertoire(Combined,group.by = 'clusters', cloneTypes=c(Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1),proportion = TRUE,
#'  label = TRUE,facet.by = NULL, na.include = FALSE,exportTable = FALSE)
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param x.axis The variable in the meta data to graph along the x.axis
#' @param label Include the number of clonotype in each category by x.axis variable
#' @param facet.by The column header used for faceting the graph
#' @param proportion Convert the stacked bars into relative proportion
#' @param na.include Visualize NA values or not.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @importFrom dplyr %>% group_by mutate
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return Stacked bar plot of counts of cells by clonotype frequency group

occupiedscRepertoire <- function(Combined,group.by=NULL, cloneTypes=c(Small = 0.001, 
                              Medium = 0.01, Large = 0.1, Hyperexpanded = 1),proportion = TRUE,
                                 label = TRUE,
                                 facet.by = NULL, 
                                 na.include = FALSE,
                                 exportTable = FALSE) {
    cloneCall <- "CTaa"
	colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

	#checkSingleObject(sc)
    meta <- grabMeta(Combined,group.by)
	Con.df <- NULL
	cell.names <- rownames(meta)
	# calculate frequency of clonotype per group.by variable
	if(proportion){
		cloneTypes <- c(None = 0, Rare = 1e-4, cloneTypes)
	}else{
		cloneTypes <- c(None = 0, Single = 1, cloneTypes)
	}
	
    data2 <- na.omit(unique(meta[,c("barcode", 'CTaa', group.by)]))
    data2 <- data2[data2[,"barcode"] %in% cell.names, ]
    data2 <- as.data.frame(data2 %>% group_by(data2[,'CTaa'], 
                data2[,group.by]) %>% summarise(Frequency = n()))
    colnames(data2)[c(1,2)] <- c('CTaa', group.by)
    x <- unique(meta[,group.by])
    for (i in seq_along(x)) {
            sub1 <- subset(meta, meta[,group.by] == x[i])
            sub2 <- subset(data2, data2[,group.by] == x[i])
            merge <- merge(sub1, sub2, by='CTaa')
            if (proportion == TRUE) {
                merge$Frequency <- merge$Frequency/length(merge$Frequency)
            }
            Con.df <- rbind.data.frame(Con.df, merge)
        } 
    nsize <- Con.df %>% group_by(Con.df[,paste0(group.by, ".x")])  %>% summarise(n = n())
	
	Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <- 
        paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], 
        ' < X <= ', cloneTypes[x], ')') }
    for (i in 2:length(cloneTypes)) { Con.df$cloneType <- 
        ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency 
        <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
    PreMeta <- unique(Con.df[,c("barcode", "CTaa", "Frequency", "cloneType")])
	dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
	
	"%!in%" <- Negate("%in%")
    PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
    rownames(PreMeta) <- PreMeta$barcode
	
	PreMeta <- PreMeta[order(match(PreMeta$barcode,meta$barcode)),]
	meta$Frequency <- PreMeta$Frequency
	meta$cloneType <- PreMeta$cloneType
    
	
    meta <- reshape2::melt(table(meta[!is.na(meta$Frequency), 
                c(group.by, facet.by, "cloneType")], useNA = "ifany"))
    if (!na.include) {
      meta <- na.omit(meta)
    }
    meta <- meta[meta$value != 0,]
    if(proportion == TRUE) {
      meta <- meta %>%
        group_by(meta[,1]) %>%
        mutate(total = sum(value), 
               prop = value/total)
      meta <- as.data.frame(meta)
    }
    if (exportTable == TRUE) {
        return(meta)
    }
    col <- length(unique(meta$cloneType))
	
    if(proportion == TRUE) {
      plot <- ggplot(meta, aes(x = factor(meta[,group.by]), y = prop, fill = cloneType)) + 
        geom_bar(stat = "identity") 
      lab <- "Proportion of Cells"
         
    } else {
      plot <- ggplot(meta, aes(x = factor(meta[,group.by]), y = value, fill = cloneType)) + 
        geom_bar(stat = "identity")+scale_fill_discrete(breaks=c('B', 'C', 'A')) 
      lab <- "Single Cells"
      
    } 
    plot <- plot + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = c(colorblind_vector(col))) + 
	  ylab(lab) + 
      theme_classic() + 
      theme(axis.title.x = element_blank())
    if (!is.null(facet.by)) {
      plot <- plot + facet_grid(.~meta[,facet.by])
    }
    if (label == TRUE) {
        plot <- plot + geom_text(aes(label = value), position = position_stack(vjust = 0.5))
      }
    return(plot)
}



#' Examine the clonal diversity of samples
#'
#' This function calculates traditional measures of diversity - Shannon, 
#' inverse Simpson, Chao1 index, abundance-based coverage estimators 
#' (ACE), and 1-Pielou's measure of species evenness by sample or group. 
#' The function automatically down samples the
#' diversity metrics using 100 boot straps The group parameter can be 
#' used to condense the individual samples. If a matrix output for 
#' the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalDiversity(combined, cloneCall = "gene")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by Variable in which to group the diversity calculation
#' @param x.axis Additional variable in which to split the x.axis
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Exports a table of the data into the global environment 
#' in addition to the visualization
#' @param n.boots number of bootstraps to downsample in order to get mean diversity
#' @importFrom stringr str_sort str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
#' @author Andrew Malone, Nick Borcherding
clonalDiversity <- function(Combined,group.by=NULL, x.axis = NULL, split.by = NULL,
                            exportTable = FALSE, n.boots = 100) {
  cloneCall <- "CTaa"
  colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

  df <- list.input.return(Combined,group.by)
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
    values <- stringr::str_sort(as.character(unique(melt[,group.by])), 
                       numeric = TRUE)
    values <- quiet(dput(values))
    melt[,group.by] <- factor(melt[,group.by], levels = values)
	
    if (x.axis == "x.axis") {
        plot <- ggplot(melt, aes(x=1, y=as.numeric(value)))
    } else {
      plot <- ggplot(melt, aes(x=melt[,x.axis], y=as.numeric(value)))
    }
    plot <- plot +
      geom_boxplot(outlier.alpha = 0) +
      geom_jitter(aes(color = melt[,group.by]), size = 3) + 
      labs(color="Group") +
      ylab("Index Score") +
      scale_color_manual(values = colorblind_vector(length(unique(melt[,group.by])))) +
    facet_wrap(~variable, scales = "free", ncol = 5) +
      theme_classic() + 
      theme(axis.title.x = element_blank())
    if (x.axis == "x.axis") { 
      plot <- plot + theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
      }
  if (exportTable == TRUE) { return(mat) }
  return(plot) 
}


#' Examining the clonal space occupied by specific clonotypes
#'
#' This function calculates the relative clonal space occupied by the 
#' clonotypes. The grouping of these clonotypes is based on the parameter 
#' split, at default, split will group the clonotypes into bins of 1:10, 
#' 11:100, 101:1001, etc. To adjust the clonotypes selected, change the 
#' numbers in the variable split. If a matrix output for the data is
#' preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param split The cutpoints for the specific clonotypes.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#'
#' @export
#' @return ggplot of the space occupied by the specific rank of clonotypes
clonalProportion <- function(Combined,group.by, split = c(10, 100, 1000, 10000, 30000, 
                        100000),exportTable = FALSE) {
    cloneCall <- "CTaa"
	colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

	Con.df <- NULL
    df <- list.input.return(Combined,group.by)
    cloneCall <- theCall(cloneCall)
    df <- checkList(df)
    df <- checkBlanks(df)
    mat <- matrix(0, length(df), length(split), dimnames = list(names(df), 
            paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))
    
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, split[-length(split)] + 1)
    for (i in seq_along(split)) {
        mat[,i] <- vapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])), 
                            FUN.VALUE = numeric(1))
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    mat_melt <- reshape2::melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", 
                        values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()
    return(plot)
}



#' @ a function to prepare input for chord diagram
#' @param sc object after combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param proportion Binary will calculate relationship unique 
#' clonotypes (proportion = FALSE) or a ratio of the group.by 
#' variable (proportion = TRUE)
#' 
#' @importFrom reshape2 dcast
#' @export
#' @return data frame of shared clonotypes between groups
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(Combined,group.by, proportion) {
    
	cloneCall <- "CTaa"
	colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

	meta <- grabMeta(Combined,group.by)
    cloneCall <- theCall(cloneCall)
    test <- meta[, c(cloneCall, group.by)]
    test <- test[!is.na(test[,cloneCall]),]
    dTest <- suppressMessages(reshape2::dcast(test, test[,cloneCall] ~ test[,group.by]))
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest2 <- dTest[-1]
    dTest2[dTest2 >= 1] <- 1
    total <- nrow(dTest)
  
    list <- list()
    for (i in seq_len(nrow(dTest2))) {
      list[[i]] <- which(dTest2[i,] > 0)
    }
    matrix_out <- matrix(ncol = ncol(dTest2), nrow = ncol(dTest2), 0)
    for (j in seq_along(list)) {
      matrix_out[list[[j]],list[[j]]] <- matrix_out[list[[j]],list[[j]]] + 1
      if (length(list[[j]]) > 1) {
        #length <- length(list[[j]])
        diag(matrix_out[list[[j]],list[[j]]]) <-  diag(matrix_out[list[[j]],list[[j]]]) - 1
      }
    }
  
    matrix_out[lower.tri(matrix_out)] <- NA

    colnames(matrix_out) <- colnames(dTest2)
    rownames(matrix_out) <- colnames(dTest2)
    
    output <- data.frame(from = rep(rownames(matrix_out), 
                        times = ncol(matrix_out)),
                        to = rep(colnames(matrix_out), each = nrow(matrix_out)),
                        value = as.vector(matrix_out),
                        stringsAsFactors = FALSE)
    output <- na.omit(output)
    
    if (proportion == TRUE) {
        output$value <- output$value/total
    } 
    return(output)
}

############### utils


