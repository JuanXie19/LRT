#' get Overall Trajectory
#'
#' @description This function infers the overall trajectory structure from \code{CombinedDataSet} using \code{\link[slingshot]}.
#'
#' @param data A class of \code{CombinedDataSet} object
#' @param reducedDim The dimensionality reduction to be used. Can be
#'   a matrix or a character identifying which element of
#'   \code{reducedDim(data)} is to be used. If multiple dimensionality
#'   reductions are present and this argument is not provided, the first element
#'   will be used by default.
#' @param start.clus (optional) character, indicates the starting cluster(s)
#'   from which lineages will be drawn.
#' @param end.clus (optional) character, indicates which cluster(s) will be
#'   forced to be leaf nodes in the graph.
#' @param dist.method (optional) character, specifies the method for calculating
#'   distances between clusters. Default is \code{"slingshot"}, see
#'   \code{\link[TrajectoryUtils]{createClusterMST}} for details.
#' @param use.median logical, whether to use the median (instead of mean) when
#'   calculating cluster centroid coordinates. Default is \code{'TRUE'}.
#' 
#' @details Given a \code{CombinedDataSet} object, this function infers the overall cell trajectory using the \code{\link[slingshot]}.
#'   To reduce noise, clones with less than 5 cells are removed before trajectory inference.
#'
#'
#' @import Seurat
#' @import dplyr
#' @import slingshot
#' @import tidyr
#' @import TrajectoryUtils
#' @import progress
#' @import SingleCellExperiment
#'
#' @examples
#' TCR <-read.csv('/PATH/TO/YOUR/scTCR-seqData/',header=T)
#' load('Mice.sub.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getOverallTrajectory(Combined,reducedDim = NULL,start.clus = NULL, end.clus = NULL, dist.method = 'slingshot', use.median = TRUE)
#'
#' @return A class of \code{TrajectoryDataSet} object
#'
#' @export 


setMethod(f = 'getOverallTrajectory',
	signature = signature('CombinedDataSet'),
	definition = function (CombinedDataSet,reducedDim = 'UMAP', start.clus = NULL, end.clus = NULL, dist.method = 'slingshot', use.median = FALSE,...){
	
	## filter the data to remove clones with less than 5 cells
	seurat <- CombinedDataSet@seurat
	temp.df <- data.frame(cell = colnames(seurat),clone=seurat$CTaa,clusters = seurat$clusters)
	temp.df1 <- temp.df %>% dplyr::group_by(clone,clusters) %>% dplyr::summarise(n=n())
	temp.df1.wide <- temp.df1%>% tidyr::pivot_wider(names_from = clusters,values_from = n,values_fill=0)
	temp.df1.wide2 <- temp.df1.wide[rowSums(temp.df1.wide[,2:ncol(temp.df1.wide)])>4,]
	
	seurat <- subset(seurat, subset = CTaa %in% temp.df1.wide2$clone)
	
	
	## infer trajectory using slingshot
	sce <- Seurat::as.SingleCellExperiment(seurat,assay='RNA')
	CLUSTERS <- SingleCellExperiment::colData(sce)$clusters
	
	sce <- slingshot(sce, reducedDim = reducedDim,start.clus = start.clus,clusterLabels = CLUSTERS,dist.method = dist.method,use.median = use.median)
	Params <-list(start.clus = start.clus, end.clus = end.clus, dist.method = dist.method, use.median = use.median)
	
		
	TrajectoryData <- new('TrajectoryDataSet',seurat = seurat,slingCurves = slingCurves(sce), slingLineages = slingLineages(sce),slingPseudotime = slingPseudotime(sce),trajectoryParams = Params)
	return(TrajectoryData)
})


