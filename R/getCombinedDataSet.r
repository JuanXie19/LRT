#' Combine scTCR-seq and scRNA-seq data 
#'
#' This function attaches the immune receptor information to metadata of the seurat object
#' @param scTCRseq The preprocessed scTCR-seq data. See details
#' @param seurat A seurat object
#'
#' @details The input scTCR-seq needs to be a data frame, and alpha chain and beta chain need to be paired.
#'  If you have multiple sets of scTCR-seq from more than one sample, the data needs to be combined into a single data frame.
#' 	It should contain a column named 'barcode' and a column named 'CTaa'. The 'barcode' column corresponds to the cell barcodes(potentially prefixed with sample name).
#'  The 'CTaa' column holds the scTCR-seq information for cells. Paired alpha-beta chain amino acid sequence is prefered, yet users are free to provide any related information that can be used to define a clonotype.
#'  For instance, you can put CDR3 alpha chain, or beta chain, or paired alpha-beta chain amino acid sequence in this column.
#'  If provide paired alpha-beta chain amino acide sequence, the alpha chain and beta chain should be separated by '_'.
#'  If you want to check the V and J gene usage pattern, you need to provide another column named 'CTgene', with alpha chain gene and beta chain gene separated by '_'. 
#'
#' @details The seurat object must contain clustering results saved as 'clusters', and must have cell barcodes saved as 'barcode'.
#'  The seurat object also needs to have an assay named 'RNA', which stores the raw count data.
#' #'  The barcodes from scTCR-seq data and seurat object must have the same naming style. For instance, if barcodes from seurat/scTCR-seq are prefixed with sample name, then the barcodes from scTCR-seq/seurat must also have prefix.
#'
#' @import Seurat
#' @import dplyr
#'
#' @export
#'
#' @return A object of class \code{CombinedDataSet} with clonotype information as well as scRNA-seq information. 
#'  The assay slot of this object is essentially an updated seurat object, with clonotype information attached.
#' 
#' @examples
#' TCR <-read.csv("/PATH/TO/YOUR/scTCR-seqData/",header=T)
#' load('Mice.sub.rda')
#' getCombinedDataSet(TCR,Mice.sub)
#'
setMethod(f='getCombinedDataSet',
		signature = signature(TCR='data.frame',seurat='Seurat'),
		definition = function(TCR,seurat,...){
			TCR <- subset(TCR, barcode %in% seurat$barcode)
			if('CTaa' %in% colnames(TCR)){
				TCR <- TCR
				}else{
				stop('No clonotype information avaialble, please provide double check and make sure it is saved in a column named as "cdr3" ')}
			
			CHECK <- identical(TCR$barcode,names(seurat$barcode))
			if(CHECK==0) stop('barcode from scTCR-seq and scRNA-seq not match, please check! ')
			TCR <- TCR[order(match(TCR$barcode,seurat$barcode)),]   # reorder the TCR dataframe to make sure we assign correct cdr3 to the seurat object
			seurat$CTaa <- TCR$CTaa
			if('CTgene' %in% colnames(TCR)){
				seurat$CTgene <- TCR$CTgene
			}
			CombinedData <-new('CombinedDataSet',seurat = seurat)
			return(CombinedData)
		}
)