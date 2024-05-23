## Script for analyzing Xbp1/Ire1KO cDC1 CITE-seq project data
## Basic barebones pipeline run on HPC for the individual samples and 
## Similar pipeline run for the initial merge objects

######################################################################
################ SPECIFY INPUT  AND OUTPUT ###########################
######################################################################

setwd('/srv/data/10x/ABC/ABC001/')
output.dir <- "results/"
samplename <- "ABC001"
experiment <- samplename
dimsToUse <- 25
resToUse <- 0.8

dir.create(output.dir)
dir.create(paste0(output.dir,"Robjects"))
dir.create(paste0(output.dir,"Plots"))
dir.create(paste0(output.dir,"Plots/RNA"))
dir.create(paste0(output.dir,"Plots/ADT"))

################################################################################
############################## LOAD PACKAGES  ##################################
################################################################################
library("tidyverse")
library("Seurat")
library("SingleCellExperiment")
library("scater")
library("scran")
library('DoubletFinder')
library('fields')
library('modes')
library("dplyr")
library("stringr")
library("readr")

runDoubletFinder <- function(object, PCs, minPCT = 1, maxPCT = 10, pN = 0.25){
  
  if(!require('DoubletFinder')){
    devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  }
  DFPredictions <- c()
  findpK <- paramSweep_v3(object, PCs = PCs, sct=T) %>%    summarizeSweep(GT=F) %>%    find.pK()
  maxScore <- findpK %>% pull('BCmetric') %>% which.max()
  pKValue <- findpK[maxScore, 'pK'] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(object@meta.data$SCT_clusters) 
  nExp_poi.min <- round((minPCT/100)*length(colnames(object)))
  nExp_poi.min.adj <- round(nExp_poi.min*(1-homotypic.prop))
  nExp_poi.max <- round((maxPCT/100)*length(colnames(object)))
  nExp_poi.max.adj <- round(nExp_poi.max*(1-homotypic.prop))
  object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue,nExp =  nExp_poi.min, reuse.pANN = F, PCs = PCs, sct = T)
  object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue,nExp =  nExp_poi.min.adj, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_", nExp_poi.min), PCs = PCs, sct = T)
  object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue,nExp =  nExp_poi.max, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_",  nExp_poi.min), PCs = PCs, sct = T)
  object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue,nExp = nExp_poi.max.adj, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_",  nExp_poi.min), PCs = PCs, sct = T)
  object@meta.data$DFPrediction <- object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nExp_poi.max)]
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.max.adj)] == 'Doublet'] <- "Low Confidence adjusted"
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.min)] == 'Doublet'] <- "High Confidence"
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.min.adj)] == 'Doublet'] <- "High Confidence adjusted"
  object@meta.data$DFPrediction <- gsub('Doublet', "Low Confidence", object@meta.data$DFPrediction)
  return(object)
}


################################################################################
################################## RNA  ########################################
################################################################################

#COUNT DATA
rawData <- Read10X_h5("rawFiles/filtered_feature_bc_matrix.h5")
rawDataRNA <- as.matrix(rawData$`Gene Expression`)

# SCE
sce<-SingleCellExperiment(list(counts=rawDataRNA))
is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
sce<- addPerCellQC(sce, subsets=list(Mito=is.mito))
metaData<-data.frame("staticNr"=colnames(rawDataRNA), 
                     "nGene"=sce$detected,
                     "nUMI"=sce$sum,
                     "percent.mito"=sce$subsets_Mito_percent, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

# FILTERING
nmad_low_feature<-5
nmad_high_feature<-5
nmad_low_UMI<-5
nmad_high_UMI<-5
nmad_high_mito<-10

feature.drop.low <- isOutlier(sce$detected, nmads=nmad_low_feature, type="lower", log=TRUE,batch = batch)
feature.drop.high <- isOutlier(sce$detected, nmads=nmad_high_feature, type="higher", log=TRUE,batch = batch)
feature.drop<-as.logical(feature.drop.low + feature.drop.high,batch = batch)
libsize.drop.low <- isOutlier(sce$sum, nmads=nmad_low_UMI, type="lower", log=TRUE,batch = batch)
libsize.drop.high <- isOutlier(sce$sum, nmads=nmad_high_UMI, type="higher", log=TRUE,batch = batch)
libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
mito.drop.high <- isOutlier(sce$subsets_Mito_percent, nmads=nmad_high_mito, type="higher",batch = batch)
mito.drop<-as.logical(mito.drop.high)

metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]

# NORMALIZE: SCRAN/SCATER
set.seed(123)
q.clust<- quickCluster(sce) 
sce <- scran::computeSumFactors(sce, cluster=q.clust)
sce <- scater::logNormCounts(sce)
saveRDS(sce, file=paste0(output.dir,"/Robjects/sce_RNA.rds"))
rawDataFiltered<-rawDataRNA[rownames(sce),colnames(sce)] 

# CREATE SEURAT OBJECT
cells.to.keep <- colnames(sce)
seuratObj_tmp <- CreateSeuratObject(rawDataFiltered, project = experiment)
seuratObj <- as.Seurat(sce, assay = "RNA", counts = "counts", data = "logcounts", project = experiment)
meta_tmp <- merge(seuratObj_tmp@meta.data, seuratObj@meta.data, by=0, all=TRUE, sort=FALSE)
rownames(meta_tmp) <- meta_tmp[,1]
meta_tmp <- meta_tmp[,-1]
seuratObj@meta.data <- meta_tmp

# SCT
seuratObj <- SCTransform(seuratObj, verbose = TRUE, new.assay.name = "SCT")

# SCALING & PCA
HVG <- VariableFeatures(seuratObj,assay="SCT")
seuratObj <- ScaleData(seuratObj, assay = "RNA")
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "SCT",
                    reduction.name = "SCT_pca",reduction.key = "sctPC_")
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "RNA",
                    reduction.name = "RNA_pca",reduction.key = "rnaPC_")

# FIND CLUSTERS
seuratObj <- FindNeighbors(object = seuratObj, reduction = "SCT_pca", dims = dimsToUse)
seuratObj <- FindClusters(object = seuratObj,  graph.name = "SCT_snn", resolution = resToUse)

# tSNE
seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "SCT", reduction = "SCT_pca", 
                       reduction.name = "SCT_tsne", reduction.key = "sctTSNE_")
# UMAP
seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "SCT", reduction ="SCT_pca",
                       reduction.name = "SCT_umap", reduction.key = "sctUMAP_")

# DOUBLET DETECTION
seuratObj <- runDoubletFinder(seuratObj, dimsToUse, minPCT = 1, maxPCT = 15, pN =0.25)

# DEG
Idents(object = seuratObj) <- 'SCT_clusters'
exp.markers <- FindAllMarkers(seuratObj,assay="SCT",slot = "data")

saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_RNA.rds"))
saveRDS(exp.markers, file = paste0(output.dir,"Robjects/ExpressionMarkers_",experiment,"_RNA.rds"))

################################################################################
################################## ADT  ########################################
################################################################################

# ADT COUNT DATA
rawDataADT <- as.matrix(rawData$`Antibody Capture`)
raw.names <- rownames(rawDataADT)
cells.use  <- intersect(colnames(rawDataADT), colnames(seuratObj))
ABS.use  <- rownames(rawDataADT)
rawDataADT <- rawDataADT[ABS.use, cells.use]
seuratObj[["ADT"]] <- CreateAssayObject(counts = rawDataADT[,colnames(seuratObj)])

# CLR transformation
seuratObj <- NormalizeData(seuratObj, assay = "ADT", normalization.method = "CLR", verbose=T)

# Informative Antibodies
seuratObj <- FindVariableFeatures(seuratObj,assay="ADT",selection.method = "vst", nfeatures=nrow(rawDataADT))
HVABs <- VariableFeatures(seuratObj, assay = "ADT")

# SCALING & PCA
seuratObj <- ScaleData(seuratObj, assay = "ADT")
seuratObj <- RunPCA(object = seuratObj, assay = "ADT", features = VariableFeatures(seuratObj, assay = "ADT"), 
                    npcs = nrow(rawDataADT), ndims.print = 1:5, nfeatures.print = 10,
                    reduction.name = "ADT_pca",reduction.key = "adtPC_")

# Find clusters
seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_pca", dims = dimsToUse)
seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = resToUse)
  
# tSNE
seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "ADT", reduction = "ADT_pca", 
                       reduction.name = "ADT_tsne", reduction.key = "adtTSNE_")

# UMAP
seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_pca",
                       reduction.name = "ADT_umap", reduction.key = "adtUMAP_")

# DEG
Idents(object = seuratObj) <- 'ADT_clusters'
exp.markers <- FindAllMarkers(seuratObj,assay="ADT")

saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_ADT.rds"))
saveRDS(exp.markers, file=paste0(output.dir,"Robjects/ExpressionMarkers_",experiment, "_ADT.rds"))

