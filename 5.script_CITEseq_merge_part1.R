## Script for analyzing Xbp1/Ire1KO cDC1 CITE-seq project data
## Part 1 of detailed pipeline run on the merge of WT and DKO cDC1 data

######################################################################
################ SPECIFY INPUT  AND OUTPUT ###########################
######################################################################
dir.create("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Final/PROCESSED_DATA/SAM_merge/")

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Final/PROCESSED_DATA/SAM_merge/")

output.dir <- "results_merge_final/"

samplename <- "SAM5-6-7-8"

experiment <- samplename

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') 

################################################################################
############################## RUN SCRIPT  ##################################### 
################################################################################

library("tidyverse")
library("Seurat")
library("SingleCellExperiment")
library("scater")
library("data.table")
library("ggpubr")
library("gridExtra")
library("fitdistrplus") 
library("dplyr")
library("plyr")
library("caret")
library("openxlsx")
library("ggplot2")
library("tibble")
library("scran")
library("limma")
library("ggrepel")
library("intrinsicDimension")
library("tidyr")
library("purrr")
library("openxlsx")
library("ggrepel")
library("clustree")
library("limma")
library("patchwork")
library("future")
library("DisneyTools")

dir.create(output.dir)
dir.create(paste0(output.dir,"Robjects"))
dir.create(paste0(output.dir,"Plots"))
dir.create(paste0(output.dir,"Plots/RNA"))

diagnostics<-list()

## Load in seuratobjects (DKO and WT subsets)
## Issue with WT subset: RNA values not log normalized correctly
## Workaround: start from WT non-subset and subset after and transfer annotations!!!
seuratObj_SAM05 <- readRDS(file="../SAM05/Robjects/seuratObj_sliced_SAM05_clint.rds")
seuratObj_SAM06 <- readRDS(file="../SAM06/Robjects/seuratObj_sliced_SAM06_clint.rds")
Idents(seuratObj_SAM06)<-seuratObj_SAM06@meta.data$annotated_clusters

## Perform merge WT with merge.data = T!!!
seuratObj_WT <- merge(seuratObj_SAM05, y = seuratObj_SAM06, 
                   add.cell.ids = c("SAM05","SAM06"), project = "SAM_WT_merge", merge.data = T)

## Load in subsets
seuratObj_WT_subset <- readRDS(file="../SAM05and06_WT/results_subset/Robjects/seuratObj_cDC1s_clean_SAM05and06_WT_2021.rds")
seuratObj_DKO_subset <- readRDS(file="../SAM07and08_DKO/results_subset/Robjects/seuratObj_cDC1s_clean_SAM07and08_DKO_2021.rds")

## Subset WT merge with colnames WT_subset object!
length(intersect(colnames(seuratObj_WT_subset),colnames(seuratObj_WT)))
length(setdiff(colnames(seuratObj_WT_subset),colnames(seuratObj_WT)))
length(setdiff(colnames(seuratObj_WT),colnames(seuratObj_WT_subset)))

seuratObj_WT_subset_norm<-subset(seuratObj_WT, cells = colnames(seuratObj_WT_subset))

## Transfer metadata and idents
Idents(seuratObj_WT_subset_norm)<-Idents(seuratObj_WT_subset)
seuratObj_WT_subset_norm@meta.data<-seuratObj_WT_subset@meta.data

## Perform merge
seuratObj <- merge(seuratObj_WT_subset_norm, y = seuratObj_DKO_subset, 
                   project = "SAM_merge", merge.data = T) 

unique(sapply(X = strsplit(colnames(seuratObj), split = "_"), FUN = "[", 1))

table(seuratObj$orig.ident)

rm(seuratObj_SAM05)
rm(seuratObj_SAM06)
rm(seuratObj_DKO_subset)
rm(seuratObj_WT)
rm(seuratObj_WT_subset)
rm(seuratObj_WT_subset_norm)
gc()

## Save idents seuratobjects
Idents(seuratObj)
seuratObj@meta.data$WT_DKO_annotation<-Idents(seuratObj)

################################
########## SCT transformation
################################

# does Normalizedata, scaledata and findvariablefeatures in one
seuratObj <- SCTransform(seuratObj, verbose = TRUE, new.assay.name = "SCT")

# SCT is now the default assay
Key(seuratObj)
names(seuratObj)

HVG <- VariableFeatures(seuratObj,assay="SCT")
length(VariableFeatures(seuratObj)) # 3000 in SCT

head(HVFInfo(seuratObj))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(output.dir,"Plots/RNA/05_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))

###########################
########## SCALING & PCA
###########################
# Scaling the RNA assay data, not the SCT!
#### to be complete we will scale the normalized data so we have an SCT alternative
seuratObj <- ScaleData(seuratObj, assay = "RNA")

# Run PCA on sct (150 PCs for more robust automatic selection)
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "SCT",
                    reduction.name = "SCT_pca",reduction.key = "sctPC_")

# Run PCA on rna normalized through scran/scater
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "RNA",
                    reduction.name = "RNA_pca",reduction.key = "rnaPC_")

names(seuratObj)
seuratObj[['SCT_pca']]@cell.embeddings[1:5,1:5]
seuratObj[['RNA_pca']]@cell.embeddings[1:5,1:5]

########################################
########## PCA PLOT
########################################
pdf(file=paste0(output.dir,"Plots/RNA/07_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
dev.off()

########################################
########## PC loadings
########################################

# plotting residual variance SCT (HVG) vs PC-loadings
# build the data-frame
loadings <- rownames_to_column(as.data.frame(seuratObj@reductions$SCT_pca@feature.loadings),var="gene")
var.load <- seuratObj@assays$SCT@meta.features %>%
  rownames_to_column(var="gene") %>%
  filter(gene %in% HVG) %>%
  # select(gene,sct.residual_variance) %>%
  left_join(loadings,by="gene")

var.load.plot <- ggplot(var.load,aes(sct.residual_variance,sctPC_2,color=sctPC_2)) +
  geom_point() +
  geom_abline(aes(slope = 0, intercept = 0)) +
  geom_label_repel(aes(label=ifelse(abs(sctPC_2)>0.1,as.character(gene),'')),color="black") +
  theme_classic() +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
var.load.plot

# Make PC plots for inspection

pc_list <- colnames(seuratObj@reductions$SCT_pca@feature.loadings)
pc_plots <- list()
for (i in 1:50) {
  labels <-ifelse(abs(var.load[,2+i])>0.1,as.character(var.load$gene),'')
  p <- ggplot(var.load,aes_string("sct.residual_variance",pc_list[[i]][1],color=pc_list[[i]][1])) +
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    geom_label_repel(aes(label=ifelse(abs(.data[[pc_list[[i]]]])>0.1,as.character(gene),'')),color="black") +
    theme_classic() +
    scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
  pc_plots[[i]] <- p
}

pdf(file=paste0(output.dir,"Plots/RNA/08a_PCloadings.pdf"))
for (i in 1:50) {
  print(pc_plots[[i]])
}
dev.off()

####################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
####################################################

# based on Pipecomp paper
int.dim <- maxLikGlobalDimEst(seuratObj@reductions$SCT_pca@cell.embeddings,
                              k=20,unbiased = TRUE, neighborhood.aggregation = 'robust')

est.PC <-round(int.dim[[1]])
est.PC
diagnostics[['est.dimsPC.maxlik']] <- est.PC

### Create PCElbowplot
png(file=paste0(output.dir,"Plots/RNA/08b_selectPC_SCT.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 50, reduction = "SCT_pca") + geom_vline(aes(xintercept=est.PC))
dev.off()

dimsToTry<-c(est.PC,seq(20,40,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "SCT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "SCT_snn", resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "SCT", reduction = "SCT_pca", 
                       reduction.name = "SCT_tsne", reduction.key = "sctTSNE_")
  tsnePlot<-DimPlot(seuratObj, reduction = "SCT_tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "SCT_tsne", label=F, group.by="ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "SCT", reduction ="SCT_pca",
                       reduction.name = "SCT_umap", reduction.key = "sctUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(30)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "SCT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "SCT_snn", resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "SCT", reduction = "SCT_pca", 
                       reduction.name = "SCT_tsne", reduction.key = "sctTSNE_")
  tsnePlot<-DimPlot(seuratObj, reduction = "SCT_tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "SCT_tsne", label=F, group.by="ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "SCT", reduction ="SCT_pca",
                       reduction.name = "SCT_umap", reduction.key = "sctUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree
Perplexity<-30
Resolution<-0.8
Perplexity_UMAP<-30

seuratObj <- FindNeighbors(object = seuratObj, reduction = "SCT_pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "SCT_snn", resolution = res)
}

pdf(file=paste0(output.dir,"Plots/RNA/10c_Clustree.pdf"))
clustree(seuratObj, prefix = "SCT_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$SCT_clusters <- seuratObj$SCT_snn_res.0.8
Idents(seuratObj) <- seuratObj$SCT_snn_res.0.8 

################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

umapPlot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, group.by= "SCT_clusters", label.size = 6)
tsnePlot<-TSNEPlot(seuratObj, reduction = "SCT_tsne", label = T, group.by= "SCT_clusters", label.size = 6)
seuratObj@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/11_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
tsnePlot
dev.off()

# Diagnostic plots
# a series of plots to assess how well the normalization has worked
# for this we will plot the library size and calculate the correlation with the libsize

# diagnostic plots
# first calculate fraction of zeros
seuratObj$fzero<- Matrix::colSums(seuratObj@assays$RNA@counts==0)/seuratObj@assays$RNA@counts@Dim[1]

diagplot <-function(object,reduction,metric){
  data <- Embeddings(object = object[[reduction]])
  data <- as.data.frame(x = data)
  assay <- "SCT"
  method <- "sctransform"
  log10UMI <- log10(FetchData(object,vars="sum"))[[1]]
  zero.fraction <- FetchData(object,vars="fzero")[[1]]
  PC_1 <- FetchData(object,vars="sctPC_1")[[1]]
  if(metric=="umi"){
    if(reduction=="SCT_pca"){
      p <- ggplot() +
        geom_point(aes(x=sctPC_1,y=sctPC_2, colour=log10UMI), data=data, size=2, shape=20)
    }
    if(reduction=="SCT_umap"){
      p <- ggplot() +
        geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=log10UMI), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
      ggtitle(paste0("UMI (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  if(metric=="zeros"){
    if(reduction=="SCT_pca"){
      p <- ggplot() +
        geom_point(aes(x=sctPC_1,y=sctPC_2, colour=zero.fraction), data=data, size=2, shape=20)
    }
    
    if(reduction=="SCT_umap"){
      p <- ggplot() +
        geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=zero.fraction), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("red","blue")) +
      ggtitle(paste0("zero fraction  (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  if(metric=="both"){
    p <- ggplot() +
      geom_point(aes(x=log10UMI,y=sctPC_1, colour=zero.fraction), data=data, size=2, shape=20)
    as <- cor(log10UMI,PC_1)
    p <- p + scale_colour_gradientn(colours = c("red","blue")) +
      ggtitle(paste0("UMI & zero fraction r=",round(as,2)," \n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  return(p +  theme_classic() )
}

pca.umi.sct <- diagplot(seuratObj,'SCT_pca','umi')
pca.zeros.sct <- diagplot(seuratObj,'SCT_pca','zeros')
pca.both.sct <- diagplot(seuratObj,'SCT_pca','both')

umap.umi.sct <- diagplot(seuratObj,'SCT_umap','umi')
umap.zeros.sct <- diagplot(seuratObj,'SCT_umap','zeros')



pdf(file=paste0(output.dir,"Plots/RNA/12a_norm_diagnostics.pdf"), width = 17*0.45, height = 12.4*0.45)
VlnPlot(seuratObj, features =c("sum")) + ggtitle("UMI count")
VlnPlot(seuratObj, features =c("fzero")) + ggtitle("fraction of zero genes")
VlnPlot(seuratObj, features =c("detected")) + ggtitle("Detected gene count")
pca.umi.sct
pca.zeros.sct
pca.both.sct
umap.umi.sct
umap.zeros.sct
dev.off()  

pdf(file=paste0(output.dir,"Plots/RNA/12b_mito_contamination.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_Mito_percent")
VlnPlot(seuratObj, features = "subsets_Mito_percent")
dev.off() 

# adding the pca.drop to the metadata of the seuratobject
seuratObj$pca.drop <- metaData$pca.drop[!metaData$final.drop]
pdf(file=paste0(output.dir,"Plots/RNA/12c_PCA_outliers.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "pca.drop")
VlnPlot(seuratObj, features = "pca.drop")
dev.off() 

pdf(file=paste0(output.dir,"Plots/RNA/13_RBC_contamination.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_RBC_percent")
VlnPlot(seuratObj, features = "subsets_RBC_percent")
dev.off()  

pdf(file=paste0(output.dir,"Plots/RNA/14_COVID_infection.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_COVID_percent")
VlnPlot(seuratObj, features = "subsets_COVID_percent")
dev.off()  

################################################################################
############################## START ADT SCRIPT  ############################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START RNA
################################################################################

dir.create(paste0(output.dir,"Plots/ADT"))

#############################
###################################################
########## CLR transformation
###################################################
# we are doing the built-in clr normalization
# clr depends on the number of cells in a sample so possibly better to do it for each sample separately
# or try for a arcsinh(5) or logicle transform

seuratObj <- NormalizeData(seuratObj, assay = "ADT", normalization.method = "CLR", verbose=T)
Key(seuratObj)
names(seuratObj)

###################################################
########## Informative Antibodies
###################################################
# we are going to calculate the HVABs to see what ABs contain the most information
# this will create an ordering of the ABs
seuratObj <- FindVariableFeatures(seuratObj,assay="ADT",selection.method = "vst", nfeatures=nrow(seuratObj@assays$ADT)) #rawDataADT
HVABs <- VariableFeatures(seuratObj, assay = "ADT")
length(VariableFeatures(seuratObj,assay = "ADT")) # should be the number of ADTs added
head(HVFInfo(seuratObj,assay = "ADT"))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj, assay = "ADT"), 10)
plot1 <- VariableFeaturePlot(object = seuratObj, assay = "ADT")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(output.dir,"Plots/ADT/01_hvg.png"), width = 850, height = 642)
plot1 + plot2
dev.off()


##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))

################################################################################
########## SCALING & PCA
################################################################################

# Scaling the ADT assay data, not the SCT!
#### to be complete we will scale the normalized data so we have an SCT alternative
seuratObj <- ScaleData(seuratObj, assay = "ADT")

# PCA only meaningfull for large panels 
seuratObj <- RunPCA(object = seuratObj, assay = "ADT", features = VariableFeatures(seuratObj, assay = "ADT"), 
                    npcs = nrow(seuratObj@assays$ADT), ndims.print = 1:5, nfeatures.print = 10, #rawDataADT
                    reduction.name = "ADT_pca",reduction.key = "adtPC_")

names(seuratObj)

########################################
########## PCA PLOT
########################################
pdf(file=paste0(output.dir,"Plots/ADT/02_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(1,3))
dev.off()

########################################
########## HEATMAP OF PCs
########################################

# plotting residual variance SCT (HVG) vs PC-loadings
# build the data-frame
loadings <- rownames_to_column(as.data.frame(seuratObj@reductions$ADT_pca@feature.loadings),var="AB")
var.load <- seuratObj@assays$ADT@meta.features %>%
  rownames_to_column(var="AB") %>%
  filter(AB %in% HVABs) %>%
  # select(gene,sct.residual_variance) %>%
  left_join(loadings,by="AB")

var.load.plot <- ggplot(var.load,aes(vst.variance.standardized,adtPC_2,color=adtPC_2)) +
  geom_point() +
  geom_abline(aes(slope = 0, intercept = 0)) +
  geom_label_repel(aes(label=ifelse(abs(adtPC_2)>0.1,as.character(AB),'')),color="black") +
  theme_classic() +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
var.load.plot

# Make PC plots for inspection
pc_list <- colnames(seuratObj@reductions$ADT_pca@feature.loadings)
pc_plots <- list()
for (i in 1:dim(seuratObj@reductions$ADT_pca@feature.loadings)[2]) {
  labels <-ifelse(abs(var.load[,2+i])>0.1,as.character(var.load$AB),'')
  p <- ggplot(var.load,aes_string("vst.variance.standardized",pc_list[[i]][1],color=pc_list[[i]][1])) +
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    geom_label_repel(aes(label=ifelse(abs(.data[[pc_list[[i]]]])>0.1,as.character(AB),'')),color="black") +
    theme_classic() +
    scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
  pc_plots[[i]] <- p
}

pdf(file=paste0(output.dir,"Plots/ADT/03a_PCloadings.pdf"))
for (i in 1:dim(seuratObj@reductions$ADT_pca@feature.loadings)[2]) {
  print(pc_plots[[i]])
}
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

# ### alternative is to use the CV algorithm of Maxime
# based on Pipecomp paper
int.dim <- maxLikGlobalDimEst(seuratObj@reductions$ADT_pca@cell.embeddings,
                              k=20,unbiased = TRUE, neighborhood.aggregation = 'robust')

est.PC <-round(int.dim[[1]])
est.PC
diagnostics[['est.dimsPC.maxlik']] <- est.PC

### Create PCElbowplot
png(file=paste0(output.dir,"Plots/ADT/03b_selectPC_SCT.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = dim(seuratObj@reductions$ADT_pca@feature.loadings)[2], reduction = "SCT_pca") + geom_vline(aes(xintercept=est.PC))
dev.off()

################################################################################
########## CLUSTER THE CELLS
################################################################################

dimsToTry<-c(est.PC,seq(10,30,by=5)) #Used to be by 2
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "ADT", reduction = "ADT_pca", 
                       reduction.name = "ADT_tsne", reduction.key = "adtTSNE_")
  tsnePlot<-DimPlot(seuratObj, reduction = "ADT_tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "ADT_tsne", label=F, group.by="ADT_snn_res.0.8", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_pca",
                       reduction.name = "ADT_umap", reduction.key = "adtUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_umap", label = F, group.by="ADT_snn_res.0.8")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

dimsToTry<-20 #est.PC
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse, assay = "ADT", reduction = "ADT_pca", 
                       reduction.name = "ADT_tsne", reduction.key = "adtTSNE_")
  tsnePlot<-DimPlot(seuratObj, reduction = "ADT_tsne", label=T, label.size = 8)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "ADT_tsne", label=F, group.by="ADT_snn_res.0.8", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_pca",
                       reduction.name = "ADT_umap", reduction.key = "adtUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_umap", label = F, group.by="ADT_snn_res.0.8")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

### Clustering: trying out clusTree
Perplexity <- 20 #est.PC
Resolution<-0.8
resolutions <- seq(0,1,by=0.1)
for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = res)
}

pdf(file=paste0(output.dir,"Plots/ADT/4c_Clustree.pdf"))
clustree(seuratObj, prefix = "ADT_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$ADT_clusters <- seuratObj$ADT_snn_res.0.8

################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, group.by= "ADT_clusters", label.size = 6)
tsnePlot<-TSNEPlot(seuratObj, reduction = "ADT_tsne", label = T, group.by= "ADT_clusters", label.size = 6)

pdf(file=paste0(output.dir,"Plots/ADT/05_tSNE_UMAP_ADT.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
tsnePlot
dev.off()

###################################
#### Cross-Modality UMAPS and tSNE
##################################
umap_rna_rna <-  DimPlot(seuratObj, reduction = "SCT_umap",group.by="SCT_clusters",label=T)  + ggtitle("RNA (SCT) clusters")
umap_rna_adt <-  DimPlot(seuratObj, reduction = "SCT_umap",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
umap_adt_adt <-  DimPlot(seuratObj, reduction = "ADT_umap",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
umap_adt_rna <-  DimPlot(seuratObj, reduction = "ADT_umap",group.by="SCT_clusters",label=T)  + ggtitle("RNA (SCT) clusters")
tsne_rna_rna <-  DimPlot(seuratObj, reduction = "SCT_tsne",group.by="SCT_clusters",label=T)  + ggtitle("RNA (SCT) clusters")
tsne_rna_adt <-  DimPlot(seuratObj, reduction = "SCT_tsne",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
tsne_adt_adt <-  DimPlot(seuratObj, reduction = "ADT_tsne",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
tsne_adt_rna <-  DimPlot(seuratObj, reduction = "ADT_tsne",group.by="SCT_clusters",label=T)  + ggtitle("RNA (SCT) clusters")

pdf(file=paste0(output.dir,"Plots/ADT/06_overview_UMAP_tSNE.pdf"), width = 17, height = 12.4)
CombinePlots(plots = list(umap_rna_rna, umap_adt_rna, umap_rna_adt, umap_adt_adt), ncol = 2)
CombinePlots(plots = list(tsne_rna_rna, tsne_adt_rna, tsne_rna_adt, tsne_adt_adt), ncol = 2)
dev.off()


#2. ADT feature plots
#################################
# plots1 <- FeaturePlot(seuratObj, features =paste0("adt_",ABS.use) ,order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
# plots1 <- lapply(X = plots1, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))

# ENKEL ADT PER 36
plots1 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[1:36]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots1 <- lapply(X = plots1, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots2 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[37:72]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots2 <- lapply(X = plots2, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots3 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[73:108]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots3 <- lapply(X = plots3, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots4 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[109:144]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots4 <- lapply(X = plots4, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots5 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[145:180]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots5 <- lapply(X = plots5, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots6 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[181:216]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots6 <- lapply(X = plots6, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots7 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[217:252]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots7 <- lapply(X = plots7, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots8 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[253:288]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots8 <- lapply(X = plots8, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))
plots9 <- FeaturePlot(seuratObj, features = paste0("adt_",ABS.use[289:length(ABS.use)]),order=T, pt.size=0.1,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots9 <- lapply(X = plots8, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))


pdf(file = paste0(output.dir,"Plots/ADT/08_Antibodies_Featureplots_q2-q98_per36.pdf"), width = 30, height = 30)
CombinePlots(plots = plots1, ncol=6)
CombinePlots(plots = plots2, ncol=6)
CombinePlots(plots = plots3, ncol=6)
CombinePlots(plots = plots4, ncol=6)
CombinePlots(plots = plots5, ncol=6)
CombinePlots(plots = plots6, ncol=6)
CombinePlots(plots = plots7, ncol=6)
CombinePlots(plots = plots8, ncol=6)
CombinePlots(plots = plots9, ncol=6)
dev.off()

####################################
#### Clustering correspondence table
####################################
clustering.table <- table(seuratObj$SCT_clusters, seuratObj$ADT_clusters)
clus.heatmap <- pheatmap::pheatmap(clustering.table, scale = "column", display_numbers = clustering.table, main="RNA (rows) by ADT (columns) clusters") 

pdf(file = paste0(output.dir,"Plots/ADT/09_ADT_RNA_heatmap.pdf"), width = 11.69, height = 8.27)
clus.heatmap
dev.off()

#################################
#### Thresholded counttables
#################################

filtered.ADT<-data.frame("AB" = rownames(seuratObj@assays$ADT@counts),
                         "ADT_counts" = rowSums(seuratObj@assays$ADT@counts)) %>%
  dplyr::filter(ADT_counts>10) 


RNA.counts <- data.frame("gene" = rownames(seuratObj@assays$RNA@counts),
                         "RNA_counts" = rowSums(seuratObj@assays$RNA@counts))

write.xlsx(filtered.ADT, file = paste0(output.dir,"Plots/ADT/10_Countsummary_ADT_",experiment, ".xlsx"))
write.xlsx(RNA.counts, file = paste0(output.dir,"Plots/ADT/10_Countsummary_RNA_",experiment, ".xlsx"))

saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_ADT_no_doublets.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, "_ADT_no_doublets.rds"))

############################
dir.create(paste0(output.dir,"Plots/HTO"))

pdf(file = paste0(output.dir,"Plots/HTO/Overview_",experiment,".pdf"), width = 10, height = 10)
DimPlot(seuratObj,group.by = "MULTI_ID", label = F, pt.size = 0.1) + ggtitle("MULTIseqDemux")
DimPlot(seuratObj,split.by = "MULTI_ID", group.by = "MULTI_ID", label = F, pt.size = 0.1, ncol = 2) + ggtitle("MULTIseqDemux")
DimPlot(seuratObj,group.by = "HTO_maxID", label = F, pt.size = 0.1) + ggtitle("MAX HTO count")
DimPlot(seuratObj,split.by = "HTO_maxID", group.by = "HTO_maxID", label = F, pt.size = 0.1, ncol = 2) + ggtitle("MAX HTO count")
dev.off()
