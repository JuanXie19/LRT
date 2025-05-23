library(scRepertoire)
library(dplyr)

dir.in <-'/path/to/scTCR-seqData/'
# setwd(dir.in)

## read raw scTCR-seq data
M1 <- read.csv(paste0(dir.in,'TCR/M1.csv'))  
M2 <- read.csv(paste0(dir.in,'TCR/M2.csv'))   
M3 <- read.csv(paste0(dir.in,'TCR/M3.csv')) 
M4 <- read.csv(paste0(dir.in,'TCR/M4.csv'))  
M5 <- read.csv(paste0(dir.in,'TCR/M5.csv'))  

## keep chains that were annotated as full-length and functional
M1 <- subset(M1, productive=='True')  
M2 <- subset(M2, productive=='True') 
M3 <- subset(M3, productive=='True')  
M4 <- subset(M4, productive=='True') 
M5 <- subset(M5, productive=='True')  

# combine the data frames into one
M_list <-list(M1,M2,M3,M4,M5)
combined <-combineTCR(M_list,ID =c('','','','',''),samples=c('M1','M2','M3','M4','M5'),cells='T-AB',removeMulti = TRUE,filterMulti = TRUE,removeNA=TRUE)  

# only keep clones with at least 2 cells
for (i in 1:5){
  combined[[i]]$barcode =gsub('__','_', combined[[i]]$barcode)
  combined[[i]] <- combined[[i]] %>%group_by(CTnt)%>%fileter(n()>1)
}

# convert list into dataframe
combined_2cells.TCR <-do.call("rbind", lapply(combined, as.data.frame)) 
write.csv(combined_2cells.TCR,'combinedTCR_2cells.csv')

############################
### scRNA-seq data #########
library(Seurat)
dir.in <-'/path/to/scRNA-seqData/'
#setwd(dir.in)

# load cell cycling genes
s.genes <- cc.genes$s.genes  
g2m.genes <- cc.genes$g2m.genes

s.genes <- paste(substring(s.genes, 1, 1), tolower(substring(s.genes, 2)), sep = "", collapse = " ") %>% strsplit(split = " ") %>% unlist()
g2m.genes <- paste(substring(g2m.genes, 1, 1), tolower(substring(g2m.genes, 2)), sep = "", collapse = " ") %>% strsplit(split = " ") %>% unlist()

# preprocess each individual dataset
PROCESS <-function(sample_num,high,count){
    data_dir <- paste0(dir.in,sample_num)
    data <- Read10X(data.dir=data_dir)
    NUM <- strsplit(sample_num,'_')%>%sapply('[',2)
    M <- CreateSeuratObject(counts = data,project=NUM)
    M[['percent.mt']] <-PercentageFeatureSet(M,pattern='^mt-')
    M <- subset(M, subset = nFeature_RNA<high & nFeature_RNA>800 ) 
    M <- subset(M,subset=nCount_RNA <count)
    M <- subset(M, subset=percent.mt <10)# 
    M <- CellCycleScoring(M, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)  
    M <- SCTransform(M,vars.to.regress=c('percent.mt','S.Score','G2M.Score','nCount_RNA')) 
    M <- RenameCells(M,add.cell.id=NUM)  # add prefix to the colnames to avoid duplicated cell names during integration
    M$cell.id <- colnames(M)
  return(M)
}

M1rna.sct <-PROCESS("GSM4815363_M1",6000,40000)
M2rna.sct <-PROCESS("GSM4815364_M2",6000,40000)
M3rna.sct <-PROCESS("GSM4815365_M3",5000,30000) 
M4rna.sct <-PROCESS("GSM4815366_M4",5000,35000)
M5rna.sct <-PROCESS("GSM4815367_M5",5000,32000)

# combine into a list
Mice.sct <- list(M1rna.sct,M2rna.sct,M3rna.sct,M4rna.sct,M5rna.sct)

remove(M1rna.sct)
remove(M2rna.sct)
remove(M3rna.sct)
remove(M4rna.sct)
remove(M5rna.sct)

## perform integration
FEATURES.sct <-SelectIntegrationFeatures(object.list=Mice.sct, nfeatures=3000)
Mice.sct <-PrepSCTIntegration(Mice.sct,anchor.features=FEATURES.sct)
anchors.sct <-FindIntegrationAnchors(object.list = Mice.sct,normalization.method='SCT',anchor.features = FEATURES.sct)  
Mouse.combined.sct <-IntegrateData(anchorset = anchors.sct,normalization.method = 'SCT')

# PCA
Mouse.combined.sct <- RunPCA(Mouse.combined.sct, verbose = FALSE)

# filter to keep cells with corresponding TCR info
TCR <-read.csv('/path/to/combinedTCR_2cells.csv',header=T)  # 

Mouse.combined.sct_sub <-subset(Mouse.combined.sct,cell.id %in%TCR$barcode)

Mouse.combined.sct_sub <- FindNeighbors( Mouse.combined.sct_sub, reduction = "pca", dims = 1:20)
Mouse.combined.sct_sub <- FindClusters( Mouse.combined.sct_sub, n.start = 10 )

Mouse.combined.sct_sub <- RunUMAP(Mouse.combined.sct_sub, reduction = "pca", dims = 1:20)
Mouse.combined.sct_sub <- RunTSNE(Mouse.combined.sct_sub,dims=1:20)

save(Mouse.combined.sct_sub,file='Mouse.combined.sct.RData')

## cluster annotation (omitted)
## subset Tcmp and Th1 clusters
Mice.sub2 <-subset(Mouse.combined.sct_sub,idents=c(0,1,5,14,16,7,8,3,2,11,9))

## modify identities
current.cluster.ids <-c(0:3,5,7:9,11,14,16)
new.cluster.ids <-c('Tcmp','Tcmp','Th1','Th1','Tcmp','Pre-Th1','Pre-Th1','Th1','Th1-inter','Tcmp','Tcmp')
names(new.cluster.ids) <-levels(Mice.sub2)
Mice.sub2 <-RenameIdents(Mice.sub2,new.cluster.ids)
Mice.sub2$clusters <-Idents(Mice.sub2)
