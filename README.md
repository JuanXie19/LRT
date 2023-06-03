# LRT
The `LRT` package is an R package for clonotype lineage inference and clustering analysis. 

Cell differentiation is a long-standing focus in developmental biology. It is widely accepted that there are no clear distinctions between cellular states, but instead a smooth trasition. Thus cells can be viewed as points along a continuum/lineage. Inference of lineage structure help understanding how cells change states. Many existing tools perform lienage inference based on single-cell transcriptomics data, owing to the fact that cells in different states express different sets of genes, and thus transriptional changes may provide clues for the dynamic biological process that cells go through. The general idea is to learn the different states that cells might go through according to the changes of gene expression, and then place cells onto the infered path. Although single-cell transcriptomics data provide detailed phenotypic information, the predicted lineage trajectories do not necessarilly reflect the true clonal relationships[1], as cells sharing same differentiation state may come from different founder cells. 

scTCR-seq data compliments the scRNA-seq data in the sense that it provides the clonotype information for each cell, and it is generally assumed that cells sharing the same clonotype share the same progenitor. By integrating scTCR-seq data with scRNA-seq data, we could have a detailed insight of transcriptional changes that occur as clonal related cells differentiate. The `LRT` package utilizes the complementary nature of these two data type  and infers cell lineage by integrative analysis of scTCR-seq data and scRNA-seq data.


## Installation

The `LRT` package can be installed directly from this repository using `devtools`.

```
devtools::install_github('JuanXie19/LRT')
```

## Vignettes

## Requirements

The `LRT` package has multiple dependencies, including `slingshot`, `tidyverse`,`Seurat` and `TrajectoryUtils`

## Usage

Below are `R` commands for construction and clustering of clonotype lineages for a pubic mouse data published by [Khatun et al](https://pubmed.ncbi.nlm.nih.gov/33201171/)

```
# load package
library(LRT)
library(slingshot)
library(tidyverse)
library(Seurat)
library(TrajectoryUtils)

# load data
TCR <-read.csv('/path/to/scTCR-seq/TCR.csv')   # scTCR-seq data
load('Mice.sub.rda')  # seurat object, with clustering labels and UMAP

# combine TCR and scRNA-seq data
Combined <-getCombinedDataSet(TCR,Mice.sub)

# use shinyClone to intuitively explore the distribution of clonotype
shinyClone(Combined)  

# infer trajectory
## it is advised to specify the starting cluster, in this example data we choose 'Tcmp'
Trajectory <- getClonotypeLineages(Combined,start.clus = 'Tcmp', end.clus = NULL, dist.method = 'simple', use.median = TRUE)

# access one of the inferred lienages
lrtLineages(Trajectory)[[1]]

# check the parameters used for lineage inference
lrtParams(Trajectory) 

# prepare input for shinyClust
G.df <- getShinyInput(Combined)

#  lineage clustering analysis
shinyClust(Trajectory,G.df)

```


## Reference
Kester, L., & van Oudenaarden, A. (2018). Single-cell transcriptomics meets lineage tracing. Cell stem cell, 23(2), 166-179.