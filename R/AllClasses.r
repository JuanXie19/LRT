#' @title Class \code{CombinedDataSet}
#'
#' @description This \code{CombinedDataSet} class holds the combined scRNA-seq and scTCR-seq information for cells, which are relevant for 
#' 	performing clonotype lineage inference and visualization.
#'
#' @slot seurat. A seurat object containing the UMAP coordinates, cluster labels and clonotype information for cells
#'  \itemize{
#'  \item{cdr3} {character. Vector of CDR3 sequence that are used to define cell clonotype }
#'  \item{cell.id}{character. Vector of cell barcodes}
#'  \item{clusters}{character or numeric. A vector of character or numeric values indicating clustering labels for cells}
#'  \item{Other terms} {other metadata stored in the Seurat object}}

#'
#' @import Seurat
#' @import dplyr
#' @export
setClass(
	'CombinedDataSet',
	slots = c(
		seurat = 'Seurat'	
	)
)

#' @title Class \code{TrajectoryDataSet}
#'
#' @description This \code{TrajectoryDataSet} class holds the inferred clonotype lineages from \code{getOverallTrajectory}.
#' @slot clonotype. A vector of clonotype that has lineages built 
#' @slot lineages. A list of lineages corresponding the the clonotypes
#' @slot trajectoryParams. Additional parameters used by \code{getOverallTrajectory}. These may specify how the minimum spanning tree on clusters was constructed:
#'  \itemize{
#'  \item{\code{start.clus}}{ character. The label of the root cluster, or a vector of cluster labels giving the root clusters of each disjoint component of th egraph}
#'  \item{\code{end.clust}}{ character. Vector of cluster labels indicating terminal clusters.}
#'  \item{\code{dist.method}}{character. Specify the method for calculating distances between clusters}
#'  \item{\code{use.median}}{logical. Whether to use the median instead of mean when calculating cluster centroid coordinates. Default is \code{'TRUE'}.}
#'  \item{Other parameters specified by \code{\link[princurve]{principal_curve}}}.}
#'
#' @import Seurat
#' @import slingshot
#' @import dplyr
#' @examples
#' # TrajectoryDataSet objects are generated from \code{getOverallTrajectory}

#' @export

setClass(
	'TrajectoryDataSet',
	slots = c(
		seurat = 'Seurat',
		slingCurves = 'list',
		slingLineages = 'list',
		slingPseudotime = 'matrix',
		trajectoryParams = 'list')		
)

