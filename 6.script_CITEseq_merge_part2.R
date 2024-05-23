## Script for analyzing Xbp1/Ire1KO cDC1 CITE-seq project data
## Part 2 of detailed pipeline run on the merge of WT and DKO cDC1 data

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('DisneyTools')
library('openxlsx')
library("ggthemes")
library('RColorBrewer')
library("UpSetR")
library('gplots')

library('tidyverse')
library('networkD3')
library('htmlwidgets')

# change the current plan to access parallelization
library(future)
plan("multiprocess", workers = 6)
plan()

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Final/PROCESSED_DATA/")

sampleName<-"SAM_merge"
sampleFolder<-paste0(sampleName,"/")

#add some subfolders
dir.create(paste0(sampleFolder,"results_merge_non_harmony"))
dir.create(paste0(sampleFolder,"results_merge_non_harmony/QC"))
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Robjects"))

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

##### Read object (original)
seuratObj <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_SAM5-6-7-8_ADT_no_doublets.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_SAM5-6-7-8_ADT_no_doublets.rds"))

##### Create new clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters

## Check individual annotation
D1<-DimPlot(seuratObj, reduction = "SCT_umap", group.by = "annotated_clusters_individual", repel = T, label = T, label.size = 2)
ggsave(D1,file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_Individual_annotation_",sampleName,".png"), width = 20,height = 12)

D2<-DimPlot(seuratObj, reduction = "SCT_umap", group.by = "annotated_clusters", repel = T, label = T, label.size = 2)
ggsave(D2,file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_Full_annotation_",sampleName,".png"), width = 20,height = 12)

D3<-DimPlot(seuratObj, reduction = "SCT_umap", group.by = "WT_DKO_annotation", repel = T, label = T, label.size = 2)
ggsave(D3,file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_Full_annotation_minim_",sampleName,".png"), width = 20,height = 12)

##############################################################################

## Save old annotation first
seuratObj@meta.data$annotated_clusters_DKOandWT<-as.factor(paste0(seuratObj@meta.data$annotated_clusters,"_separate"))
seuratObj@meta.data$sliced_clusters_DKOandWT<-as.factor(paste0(seuratObj@meta.data$sliced_clusters,"_separate"))

# 5 distinct groups on SCT: res 0.8
# - Cl12: pre-cDC1s
# - Cl11: end of prolif
# - Cl10+16: IFN1+ Res cDC1s
# - Cl(9+13)+17: Mig cDC1s
# - Cl(14+8)+6: Prolif cDC1s
# - Rest of clusters: Res cDC1s
# - Cl((2+7)+0)+1: Earl imm
# - Cl(3+5)+(4+15): Late imm


# Check SCT 1.0 clusters!! 
DimPlot(seuratObj, reduction = "SCT_umap", group.by = "SCT_snn_res.1", label = T, label.size = 4)

## Check ADT:
diagplot2 <-function(object,reduction,metric){
  data <- Embeddings(object = object[[reduction]])
  data <- as.data.frame(x = data)
  assay <- "SCT"
  method <- "sctransform"
  log10ADT <- log10(FetchData(object,vars="nCount_ADT"))[[1]]
  FeatureADT <- FetchData(object,vars="nFeature_ADT")[[1]]
  if(metric=="ADT1"){
    if(reduction=="SCT_umap"){
      p <- ggplot() +
        geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=log10ADT), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
      ggtitle(paste0("log10 ADT count (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  if(metric=="ADT2"){
    if(reduction=="SCT_umap"){
      p <- ggplot() +
        geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=FeatureADT), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
      ggtitle(paste0("ADT feature  (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  return(p +  theme_classic() )
}

umap.ADT1.sct <- diagplot2(seuratObj,'SCT_umap','ADT1')
umap.ADT2.sct <- diagplot2(seuratObj,'SCT_umap','ADT2')

ggsave(umap.ADT1.sct, file=paste0(sampleFolder,"results_merge_non_harmony/QC/ADT_count_log10_UMAP_",sampleName,".png"), height = 7, width = 10, dpi = "retina")
ggsave(umap.ADT2.sct, file=paste0(sampleFolder,"results_merge_non_harmony/QC/ADT_feature_count_UMAP_",sampleName,".png"), height = 7, width = 10, dpi = "retina")

## Save new clustering
seuratObj@meta.data$sliced_clusters <- seuratObj@active.ident #Sliced clustering
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident #Sliced clustering

##### Create annotated clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
# logTable<-as.matrix(seuratObj[['RNA']]@data)
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

##### Read sliced object
seuratObj <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_sliced_",sampleName,"_clint.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_sliced_",sampleName,"_clint.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

################################################################################
########## PLOTS
################################################################################
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

## New function
drawUMI_mitoPlot_new<-function(coordsTable, reductionType, clusterMatrix, columnName, titleInfo){
  
  columnNr<-which(colnames(clusterMatrix)==columnName)
  
  p <- ggplot()+
    geom_point(aes(x=sctTSNE_1,y=sctTSNE_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  }
  
  p<-p +
    scale_colour_gradientn(colours = c("darkblue","cyan","green","yellow","orange","darkred")) +
    ggtitle(paste0(titleInfo," (",reductionType,")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  return(p)
}

########## UMI plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'nCount_RNA',"UMI")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"UMI")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge_non_harmony/QC/11a_UMI.png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_merge_non_harmony/QC/11b_percMito.png"), width = 20)

########## PCA plot ##########
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,3))
dev.off()

# DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "seurat_clusters")
U1<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on SCT UMAP")
U2<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on ADT UMAP")

#ADT clustering on SCT UMAP and vice versa
U3<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on SCT UMAP")
U4<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on ADT UMAP")


T1<-DimPlot(seuratObj, reduction = "SCT_tsne", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on SCT tSNE")
T2<-DimPlot(seuratObj, reduction = "ADT_tsne", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on ADT tSNE")

#ADT clustering on SCT tSNE and vice versa
T3<-DimPlot(seuratObj, reduction = "SCT_tsne", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on SCT tSNE")
T4<-DimPlot(seuratObj, reduction = "ADT_tsne", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on ADT tSNE")

################################################################################
########## ANNOTATION
################################################################################
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Annotation/"))

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/1_annotation_",sampleName,".pdf"), width = 15)
grid.arrange(U1, U3, ncol=2)
grid.arrange(U2, U4, ncol=2)
grid.arrange(T1, T3, ncol=2)
grid.arrange(T2, T4, ncol=2)
dev.off()


################################################################################
########## COLOR CELLS ACCORDING TO ADT on SCT + vice versa
################################################################################
DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="SCT_clusters")
DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="ADT_clusters")

#ADT clusters
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/1_annotation_Color_ADT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$ADT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$ADT_clusters==i),])))
  C1<-C1+ggtitle(paste0("ADT_cluster_",i))
  print(C1)
}
dev.off()

#SCT clusters
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/1_annotation_Color_SCT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$SCT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$SCT_clusters==i),])))
  C1<-C1+ggtitle(paste0("SCT_cluster_",i))
  print(C1)
}
dev.off()

################################################################################
########## CHECK DE GENES 
################################################################################
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Feature_plots"))

## Markers article Cédric ##
# pDCs: Siglech, Tcf4, Ccr9, Bst2
# cDCs: Flt3, Dpp, Zbtb46
# cDCs were further subdivided into:
# - ‘‘non-migratory’’ (Ccr7 lo Nudt17 lo ) cDC1s (Cd81 + Gcsam + Xcr1 + Irf8 + ) and cDC2s (Sirpa + S100a4 + CD209a + ) 
# - ‘‘migratory’’ (Ccr7 hi Nudt17 hi profile) cDC1s (Cd81 + Gcsam + Laptmb4 + Ncoa7 + ; 
#   the latter two genes were chosen because of downregulation of Xcr1 mRNA).
#   and cDC2s S100a4 + Anxa3 + Cdc42ep3 + ). 
# - A small group of cDCs highly expressed genes associated with cell proliferation (Mki67 + Stmn1 + ), 
#   confirming reports that cDCs can proliferate in the lungs. 
# - MCs (expression of Mafb, Mertk, and Fcgr1 and lack of expression of Flt3, Zbtb46, and Dpp4). 
# - In the cDC2 compartment, inf-cDC2s were identified based on expression of Fcgr1 and other genes (e.g., Irf8 and Irf4). 
# The top DEG of the inf-cDC2s subset compared with non-migratory cDC2s and migratory cDC2s, 
# the most discriminative genes in inf- cDC2s were related to a type I IFN signature 
# (e.g., Stat1, Iigp1, Rsad2, ifi205, Isg20, Ifit1, and Ifit3), 
# genes encoding for pro-in- flammatory cytokines (e.g., Ccl5, Cxcl9, and Cxcl10), 
# and genes promoting Th1 cell differentiation (e.g., Il12b), 
# whereas genes involved in Th2 cell development appeared to be downregulated (e.g., Irf4, Klf4, and Mgl2).

# CD209 mRNA was only found on cDC2s. Cell surface CD88 and CX3CR1 are potentially much better and stable positive discriminators of MCs, 
# particularly in combination with CD26 to identify cDCs!! We no longer recommend exclusive use of CD64 or MAR-1 to separate cDCs from MCs

# In mock-infected mice, Xcr1, Clec9a, Cadm1, and Irf8 were highly specific for cDC1s and Mertk, Emr1 (F4/80), and Fcgr1 (CD64) for MCs. 
# No cDC2-specific genes could be defined because most genes were shared with cDC1s (Dpp4 [CD26] and Zbtb46) 
# or MCs (Sirpa, Itgam, Csfr1, and Clec4a3). 
# Upon infection, a small fraction of the genes encoding for (surface) markers remained 
# cDC1 specific (e.g., Xcr1, Cadm1, Clec9a, Gcsam, Zdhhc2, Tlr3, and Itgae) or MC specific (e.g., Mertk, Mafb, and Emr1). 
# As for cDC2s in the steady state, there were no genes specific to inf-cDC2s. 
# All genes that were upregulated in inf-cDC2s compared with steady-state cDC2s were 
# shared with cDC1s (e.g., Irf8, Il12b, Il15, Cd80, Cd86, and Tnfrsf4), MCs (e.g., Fcgr1 and Tlr7), or both subsets.

##### cDC marker: Flt3, Dpp4, Zbtb46
FeaturePlot(object = seuratObj, features = "Flt3", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Dpp4", cols = c("grey", "blue"),
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Zbtb46", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

Features<-c("Flt3", "Dpp4", "Zbtb46")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"cDCS"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

##### Resident markers: (Ccr7 lo Nudt17 lo ) cDC1s (Cd81 + Gcsam + Xcr1 + Irf8 + ) 
FeaturePlot(object = seuratObj, features = "Ccr7", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Nudt17", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Cd81", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Gcsam", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Xcr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Irf8", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

Features<-c("Ccr7", "Nudt17","Cd81", "Gcsam","Xcr1", "Irf8")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Res_cDC1s"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


##### Resident markers: (Ccr7 lo Nudt17 lo ) cDC2s (Sirpa + S100a4 + CD209a + )
FeaturePlot(object = seuratObj, features = "Ccr7", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Nudt17", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Sirpa", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "S100a4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Cd209a", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

Features<-c("Ccr7", "Nudt17","Sirpa", "S100a4","Cd209a")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Res_cDC2s"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


# Migratory markers: (Ccr7 hi Nudt17 hi profile) cDC1s (Cd81 + Gcsam + Laptm4b + Ncoa7 +)
FeaturePlot(object = seuratObj, features = "Ccr7", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Nudt17", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Cd81", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Gcsam", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Laptm4b", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Ncoa7", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

Features<-c("Ccr7", "Nudt17","Cd81", "Gcsam","Laptm4b","Ncoa7")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mig_cDC1s"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

# Migratory markers: (Ccr7 hi Nudt17 hi profile) cDC2s S100a4 + Anxa3 + Cdc42ep3 )
FeaturePlot(object = seuratObj, features = "Ccr7", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Nudt17", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Anxa3", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "S100a4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
FeaturePlot(object = seuratObj, features = "Cdc42ep3", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

Features<-c("Ccr7", "Nudt17","Anxa3", "S100a4","Cdc42ep3")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Mig_cDC2s"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

# Proliferating DCs (Mki67 + Stmn1 + )
FeaturePlot(object = seuratObj, features = "Mki67", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Stmn1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Birc5", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 

Features<-c("Mki67", "Stmn1","Birc5")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Proliferating_Cells"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

# inf-cDC2s  Fcgr1 and other genes (e.g., Irf8 and Irf4). -> inflammatory model (Steady State ????????)
# type I IFN signature (Stat1, Iigp1, Rsad2, ifi205, Isg20, Ifit1, and Ifit3), 
# pro-in- flammatory cytokines (e.g., Ccl5, Cxcl9, and Cxcl10), 
# Th1 cell differentiation (e.g., Il12b), 
# Th2 cell development downregulated (e.g., Irf4, Klf4, and Mgl2)
FeaturePlot(object = seuratObj, features = "Fcgr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Irf8", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Irf4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative

Features<-c("Fcgr1", "Irf8","Irf4")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"InfCDC2s_general"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


FeaturePlot(object = seuratObj, features = "Ccl5", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Cxcl9", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Cxcl10", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 

Features<-c("Ccl5", "Cxcl9","Cxcl10")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"InfCDC2s_proInflCyt"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


FeaturePlot(object = seuratObj, features = "Il12b", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Irf4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Klf4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative
FeaturePlot(object = seuratObj, features = "Mgl2", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Negative

Features<-c("Il12b", "Irf4","Klf4","Mgl2")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"InfCDC2s_Th1_Th2"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


# CD209 mRNA only found on cDC2s. CD88 and CX3CR1 discriminators of MCs, in combination with CD26 to identify cDCs!! 
# Extra: CD64 or MAR-1 to separate cDCs from MCs
FeaturePlot(object = seuratObj, features = "Cd209a", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "C5ar1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Cd88
FeaturePlot(object = seuratObj, features = "Cx3cr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Dpp4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Cd26
FeaturePlot(object = seuratObj, features = "Fcgr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Cd64
# FeaturePlot(object = seuratObj, features = "Mar1", cols = c("grey", "blue"), 
#             reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #L2a???

# mock-infected mice, Xcr1, Clec9a, Cadm1, and Irf8 for cDC1s and Mertk, Emr1 (F4/80), and Fcgr1 (CD64) for MCs. 
# No cDC2-specific genes because most genes were shared with cDC1s (Dpp4 [CD26] and Zbtb46) or MCs (Sirpa, Itgam, Csfr1, and Clec4a3).
FeaturePlot(object = seuratObj, features = "Xcr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Clec9a", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Cadm1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Irf8", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 

Features<-c("Clec9a", "Cadm1","Xcr1","Irf8")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Extra_SScDC1s_general"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


FeaturePlot(object = seuratObj, features = "Mertk", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Adgre1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Emr1 or F4/80
FeaturePlot(object = seuratObj, features = "Fcgr1", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Cd64

FeaturePlot(object = seuratObj, features = "Dpp4", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Zbtb46", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Sirpa", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Itgam", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 
FeaturePlot(object = seuratObj, features = "Csf1r", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) #Csfr1
FeaturePlot(object = seuratObj, features = "Clec4a3", cols = c("grey", "blue"), 
            reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5) 

Features<-c("Dpp4", "Zbtb46","Sirpa","Itgam","Csf1r","Clec4a3")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Extra_SScDC2s_NonSpecific"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_tsne", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_tSNE_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

####### Extra: pre-DCs
Features<-c("Ly6c2","Ly6C","Siglech","Siglec","Irf8", "Irf4")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"pre-DCs1"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Tcf4","Bst2","Ccr2","Cd209a","S100a4","Sirpa")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"pre-DCs2"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")


####### Extra: annotation markers
Features<-c("Egr1","Egr2","Egr3","Nr4a2","Nr4a3")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Annot_markers1"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Creld2","Pdia6","Hspa5","Pdia4","Hspa1a","Hspa1b")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Annot_markers2"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("Ccr2","S100a10","Ccl3","Ccl4","Cxcl9","Cxcl10")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Annot_markers3"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

Features<-c("CD207-mh","Il12b","Ccr7","Cd63","ESAM","IA-IE")
Dimensions<-(length(Features)/2)*5
Cellpop_name<-"Annot_markers4"
F1<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1)
F2<-FeaturePlot(object = seuratObj, features = Features, cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', order = T, pt.size = 1)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")
ggsave(F2, file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Feature_plot_",Cellpop_name,"_UMAP_ordered_",sampleName,".png"), height = Dimensions, width = Dimensions+5, dpi = "retina")

#########################################################################################################################################
################################################################################
########## CHECK PERFORMANCE ADT MARKERS (DC specific + all)
################################################################################
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Test_ADT_expression"))

## Load in whitelist
ADT_SAM_df<-read.xlsx("../../CITEseq_Test/PROCESSED_DATA/PanelSAM_Clint.xlsx")

## Overview of genes used for annotation
wanted.genes<-unique(c("Siglech", "Tcf4","Ccr9", "Bst2","Flt3", "Dpp4", "Zbtb46", "Ccr7", "Nudt17","Cd81", "Gcsam","Xcr1", "Irf8", "Ccr7", 
                       "Nudt17","Sirpa", "S100a4","Cd209a", "Ccr7", "Nudt17","Cd81", "Gcsam","Laptm4b","Ncoa7",
                       "Ccr7", "Nudt17","Anxa3", "S100a4","Cdc42ep3", "Mki67", "Stmn1","Birc5", "Mafb", "Mertk",
                       "Fcgr1","Flt3", "Zbtb46","Dpp4", "Fcgr1", "Irf8","Irf4", "Stat1", "Iigp1", "Rsad2", "Ifi205", 
                       "Isg20", "Ifit1", "Ifit3", "Ccl5", "Cxcl9","Cxcl10", "Il12b", "Irf4","Klf4","Mgl2", "Clec9a", 
                       "Cadm1","Xcr1","Irf8", "Dpp4","Zbtb46", "Sirpa", "Itgam", "Csf1r", "Clec4a3",
                       "Adgre1", "Csf1r","Fcgr1","Klrb1c", "Gzmb", "Plac8", "Ly6c2", "Ccr2",
                       "Esam","Itgae","Cd207","Sell","H2-Eb1")) #Added extra
wanted.features<-c()

## Check which of those genes have matching ADT markers
for (i in wanted.genes) {
  gene_wanted<-paste0("\\<",i,"\\>")
  if(length(grep(gene_wanted,ADT_SAM_df$gene))==1) {
    wanted.features<-c(wanted.features,i,ADT_SAM_df[which(ADT_SAM_df$gene==i),2])
  }
}

## Combine in one plot: protein (ADT) levels are on the right, and RNA levels are on the left
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Test_ADT_expression/Test_ADT_expression_important_markers_",sampleName,".pdf"), width = 20)
for (i in 1:(length(wanted.features)/2)) {
  F1<-FeaturePlot(seuratObj, features = wanted.features[c(2*i-1,2*i)], min.cutoff = "q2", max.cutoff = "q98", ncol = 2)
  print(F1)
}
dev.off()


## Overview of all ADT markers with gene linked in whitelist
ADT_SAM_df<-ADT_SAM_df[order(ADT_SAM_df$gene),]
wanted.genes.full<-ADT_SAM_df$gene[!is.na(ADT_SAM_df$gene)]
wanted.ADT.full<-ADT_SAM_df$whitelist_name[!is.na(ADT_SAM_df$gene)]
wanted.features.full<-c()

## Get the matching ADT and genes in a vector
for (i in 1:length(wanted.genes.full)) {
  wanted.features.full<-c(wanted.features.full,wanted.genes.full[i],wanted.ADT.full[i])
}

# Combine all markers in one plot: protein (ADT) levels are on the right, and RNA levels are on the left 
## (split up in multiple files for easier loading)
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Test_ADT_expression/Test_ADT_expression_full_",sampleName,".pdf"), width = 20)
for (i in 1:(length(wanted.features.full)/2)) {
  F1<-FeaturePlot(seuratObj, features = wanted.features.full[c(2*i-1,2*i)], min.cutoff = "q2", max.cutoff = "q98", ncol = 2)
  print(F1)
}
dev.off()

#########################################################################################################################################
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Marker_lists"))

Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

### Find ADTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_ADTclus$cluster)
saveRDS(ADTMarkers_ADTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerADTcluster']]<-paste0(table(ADTMarkers_ADTclus$cluster)," ADT markers for ADT cluster ",rownames(table(ADTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_ADTclus<-max(as.numeric(names(table(ADTMarkers_ADTclus$cluster))))
totalNrADTclusters_ADTclusPlusOne<-totalNrADTclusters_ADTclus+1
ADTmarkersList_ADTclus<-list()

for(i in 1:totalNrADTclusters_ADTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-ADTMarkers_ADTclus[ADTMarkers_ADTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_ADTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList_ADTclus)<-paste0("ADTcluster",0:totalNrADTclusters_ADTclus)

### Write to Excel
write.xlsx(ADTmarkersList_ADTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_ADTclus$cluster)
saveRDS(SCTMarkers_ADTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerADTcluster']]<-paste0(table(SCTMarkers_ADTclus$cluster)," SCT markers for ADT cluster ",rownames(table(SCTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_ADTclus<-max(as.numeric(names(table(SCTMarkers_ADTclus$cluster))))
totalNrSCTclusters_ADTclusPlusOne<-totalNrSCTclusters_ADTclus+1
SCTmarkersList_ADTclus<-list()

for(i in 1:totalNrSCTclusters_ADTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-SCTMarkers_ADTclus[SCTMarkers_ADTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_ADTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(SCTmarkersList_ADTclus)<-paste0("ADTcluster",0:totalNrSCTclusters_ADTclus)

### Write to Excel
write.xlsx(SCTmarkersList_ADTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters

ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-max(as.numeric(names(table(ADTMarkers_SCTclus$cluster))))
totalNrADTclusters_SCTclusPlusOne<-totalNrADTclusters_SCTclus+1
ADTmarkersList_SCTclus<-list()

for(i in 1:totalNrADTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrADTclusters_SCTclus)

### Write to Excel
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-max(as.numeric(names(table(SCTMarkers_SCTclus$cluster))))
totalNrSCTclusters_SCTclusPlusOne<-totalNrSCTclusters_SCTclus+1
SCTmarkersList_SCTclus<-list()

for(i in 1:totalNrSCTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(SCTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrSCTclusters_SCTclus)

### Write to Excel
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################

U1 <- DimPlot(seuratObj, reduction = "SCT_umap", label = T, group.by = "sliced_clusters", label.size = 4)
ggsave(U1, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_UMAP_sliced_",sampleName,".png"), height = 15, width = 20, dpi = "retina")

########################################
##### Markers annotated clusters
########################################
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident 
levels(seuratObj@meta.data$annotated_clusters) <- c("Early imm Res cDC1s 1","Hspa1+CD207+ Res cDC1s","Early imm Res cDC1s 2",
                                                    "Nr4a3+CD207+ Res cDC1s","Egr1+Egr3+ Res cDC1s",
                                                    "Ccl3+Ccl4+ Res cDC1s","Prolif cDC1s 1","Apoe+ Res cDC1s","Prolif cDC1s 2",
                                                    "Mig cDC1s","Cxcl9+Cxcl10+ Res cDC1s","Prolif cDC1s 4","pre-cDC1s","pre-Mig cDC1s",
                                                    "Prolif cDC1s 3","Lower quality Res cDC1s",'IFN1+ Res cDC1s',"Mig cDC1/pre-cDC2 doublets")

U2 <- DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U2, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_UMAP_annotated_",sampleName,".png"), height = 15, width = 20, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters
colorSomeCells(clusterMatrix,umapTable, WhichCells(seuratObj, idents=c("Lower quality Res cDC1s","Mig cDC1/pre-cDC2 doublets")))

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterannotated']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){

  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

######################################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTclusterannotated']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-names(table(ADTMarkers_SCTclus$cluster))
ADTmarkersList_SCTclus<-list()

for(i in totalNrADTclusters_SCTclus){

  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

######################################################
### Make heatmap for annotated clusters
SCTMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObj.small <- subset(seuratObj, downsample = 300)

## Heatmap SCT markers on SCT clusters
top10 <- SCTMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results_merge_non_harmony/Heatmaps/Heatmap_SCTmarkersList_Annotatedclus_",sampleName,".png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_merge_non_harmony/Heatmaps/Heatmap_SCTmarkersList_Annotatedclus_",sampleName,".pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
dev.off()

################################################################################################################
##Extra: split by treatment
U_split<-DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by = "orig.ident", label.size = 4)
ggsave(U_split, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/2_UMAP_split_by_sample_",sampleName,".png"), height = 15, width = 20, dpi = "retina")

# create a dataset
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@active.ident
Aggr <- rep(sampleName,length(cluster)) 

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

barplotAggr(seuratObj, listLabels)

# Stacked
png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/3_SampleDistribution_ggplot2_annot_2_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/3_SampleDistribution_ggplot2_annot_3_",sampleName,".png"), width = 2000, height = 1500, res = 300)
ggplot(data2, aes(fill=cluster, y=Freq, x=Aggr)) + theme_bw() + 
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="SAM05"),]
data_split2<-data[which(data$Sample=="SAM06"),]
data_split3<-data[which(data$Sample=="SAM07"),]
data_split4<-data[which(data$Sample=="SAM08"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4)

png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/3_SampleDistribution_ggplot2_annot_1_",sampleName,"_adjusted.png"), width = 2000, height = 1500, res = 300)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

#############################
### Extra detail clusters
############################
Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters

detail1_15vs4<-FindMarkers(seuratObj, ident.1 = 15, ident.2 = 4, min.pct = 0.10,
                            min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

detail2_17vs9<-FindMarkers(seuratObj, ident.1 = 17, ident.2 = 9, min.pct = 0.10,
                           min.diff.pct=0.25, logfc.threshold = 0.30, only.pos = FALSE)

##### Create list
listDEgenesExtra<-tibble::lst(detail1_15vs4, detail2_17vs9)

##Add geneSymbol in column (for the export)
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                             x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
##Sort on logFC
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/detailClusters_sliced_",sampleName,".rds"))

##write to Excel
write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results_merge_non_harmony/Marker_lists/detailClusters_sliced_",sampleName,".xlsx"))


## Revert idents to annotation
Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

############################################################################################################################
# # Extra: subset seuratobj!!
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4)

# Remove the bad clusters!! (doublets + Lower quality)
# Cluster 11 ADT = Cluster17 RNA
# Doublets with pre-cDC2s?? Siglech and Ly6c2 expr??
# Not clustered separately in previous UMAPs due to not enough cells??

seuratObjNew<-subset(seuratObj, idents =c(levels(Idents(seuratObj))[c(-16,-18)]))

## Create UMAP of subset!!
U3 <- DimPlot(seuratObjNew, reduction = "SCT_umap", group.by = "annotated_clusters", label = T, repel = T, label.size = 4)
ggsave(U3, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/4_Annotated_UMAP_subset_",sampleName,".png"), width = 20, height= 15, dpi = "retina")

## Create minimal annotation of subset (Old minimal annotation in WT_DKO_annotation!!!)
seuratObjNew@meta.data$annotated_clusters_minim <- seuratObjNew@meta.data$annotated_clusters
seuratObjNew@meta.data$annotated_clusters_minim <-as.factor(as.character(seuratObjNew@meta.data$annotated_clusters_minim)) # Redo factor!!
levels(seuratObjNew@meta.data$annotated_clusters_minim) <- c("Late immature cDC1s","Late immature cDC1s","Early mature cDC1s",
                                                             rep("Early immature cDC1s",2),"Late immature cDC1s","Late immature cDC1s","IFN+ Early mature cDC1s",
                                                             'Late mature cDC1s',"Late immature cDC1s","pre-cDC1s","Early mature cDC1s",rep("Proliferating cDC1s",4))

seuratObjNew@meta.data$annotated_clusters_minim<-factor(seuratObjNew@meta.data$annotated_clusters_minim,
                                                        as.character(c("pre-cDC1s","Proliferating cDC1s","Early immature cDC1s","Late immature cDC1s",
                                                                       "IFN+ Early mature cDC1s","Early mature cDC1s",'Late mature cDC1s'))) #reorder levels

Colorset<-c(brewer.pal(10,"Paired"),"magenta","Darkred","turquoise2")

U4 <- DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters_minim", label.size = 4, cols = Colorset)
ggsave(U4, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/4_UMAP_annotated_subset",sampleName,"_minimal_annotation.png"), height = 15, width = 20, dpi = "retina")

## Use minimal annotation instead of detaile dannotation (diff with WP vs RP diff markers!!!)
Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters_minim

### Create new clusters: split on source
seuratObjNew@meta.data$newClustersTmp<-Idents(seuratObjNew)
seuratObjNew@meta.data$Genotype<-as.character(seuratObjNew@meta.data$orig.ident)
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Genotype == "SAM05"),"Genotype"]<-"WT"
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Genotype == "SAM06"),"Genotype"]<-"WT"
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Genotype == "SAM07"),"Genotype"]<-"DKO"
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Genotype == "SAM08"),"Genotype"]<-"DKO"
seuratObjNew@meta.data$newClusters<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$Genotype)
head(seuratObjNew@meta.data)

### Use the new clusters
seuratObjNew@meta.data$newClusters<- as.factor(seuratObjNew@meta.data$newClusters) #reorder levels
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters

DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, label.size = 4)

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/4_Annotated_UMAP_split_",sampleName,"_2.pdf"), width = 20, height = 12)
DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, label.size = 4)
dev.off()

########################################################################################################################


#################################################################
########## GET DIFF MARKERS DKO vs WT ##########
#################################################################

getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2,
                             min.pct = 0.10) # No min diff pct!! , min.diff.pct = 0.15 or logFC logfc.threshold = 0.30,
  markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
  
  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001
  
  markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
                     markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}

# ########## 2. GET MARKERS (everything!!) ##########
# getDEgenes<-function(ident1, ident2){
#   markersDiff <- FindMarkers(seuratObjNew, ident.1 = ident1, ident.2 = ident2,
#                              logfc.threshold = 0.01, min.pct = 0.01) #0.30, 0.10, min.diff.pct = 0.15
#   #markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
#   markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
# 
#   markersDiff$geneSymbol<-rownames(markersDiff)
#   markersDiff$pct.1<-markersDiff$pct.1+0.001
#   markersDiff$pct.2<-markersDiff$pct.2+0.001
# 
#   markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
#                      markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
#   markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
#   return(markersDiff)
# }

#### Get diff markers WP vs RP #####
levels(Idents(seuratObjNew))

precDC1s_DKOvsWT<-getDEgenes("pre-cDC1s_DKO","pre-cDC1s_WT")
precDC1s_DKOvsWT<-precDC1s_DKOvsWT[order(precDC1s_DKOvsWT$avg_logFC,decreasing = T),]
head(precDC1s_DKOvsWT)
dim(precDC1s_DKOvsWT)

Proliferating_cDC1s_DKOvsWT<-getDEgenes("Proliferating cDC1s_DKO","Proliferating cDC1s_WT")
Proliferating_cDC1s_DKOvsWT<-Proliferating_cDC1s_DKOvsWT[order(Proliferating_cDC1s_DKOvsWT$avg_logFC,decreasing = T),]
head(Proliferating_cDC1s_DKOvsWT)
dim(Proliferating_cDC1s_DKOvsWT)

Early_imm_DKOvsWT<-getDEgenes("Early immature cDC1s_DKO","Early immature cDC1s_WT")
Early_imm_DKOvsWT<-Early_imm_DKOvsWT[order(Early_imm_DKOvsWT$avg_logFC,decreasing = T),]
head(Early_imm_DKOvsWT)
dim(Early_imm_DKOvsWT)

Late_imm_DKOvsWT<-getDEgenes("Late immature cDC1s_DKO","Late immature cDC1s_WT")
Late_imm_DKOvsWT<-Late_imm_DKOvsWT[order(Late_imm_DKOvsWT$avg_logFC,decreasing = T),]
head(Late_imm_DKOvsWT)
dim(Late_imm_DKOvsWT)

IFN1_early_mat_DKOvsWT<-getDEgenes("IFN+ Early mature cDC1s_DKO","IFN+ Early mature cDC1s_WT")
IFN1_early_mat_DKOvsWT<-IFN1_early_mat_DKOvsWT[order(IFN1_early_mat_DKOvsWT$avg_logFC,decreasing = T),]
head(IFN1_early_mat_DKOvsWT)
dim(IFN1_early_mat_DKOvsWT)

Early_mat_DKOvsWT<-getDEgenes("Early mature cDC1s_DKO","Early mature cDC1s_WT")
Early_mat_DKOvsWT<-Early_mat_DKOvsWT[order(Early_mat_DKOvsWT$avg_logFC,decreasing = T),]
head(Early_mat_DKOvsWT)
dim(Early_mat_DKOvsWT)

Late_mat_DKOvsWT<-getDEgenes("Late mature cDC1s_DKO","Late mature cDC1s_WT")
Late_mat_DKOvsWT<-Late_mat_DKOvsWT[order(Late_mat_DKOvsWT$avg_logFC,decreasing = T),]
head(Late_mat_DKOvsWT)
dim(Late_mat_DKOvsWT)

##add to list
listDiffMarkers<-tibble::lst(precDC1s_DKOvsWT,Proliferating_cDC1s_DKOvsWT,Early_imm_DKOvsWT,Late_imm_DKOvsWT,
                             IFN1_early_mat_DKOvsWT,Early_mat_DKOvsWT,Late_mat_DKOvsWT)

lapply(listDiffMarkers, dim)
listDiffMarkers<-lapply(listDiffMarkers,function(x){x<-x[order(x$score, decreasing=T),]})

#Check settings
saveRDS(listDiffMarkers,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/markersDiffSamples_Full_",sampleName,".rds"))

### Write to Excel
write.xlsx(listDiffMarkers, file = paste0(sampleFolder,"results_merge_non_harmony/Marker_lists/summaryDiffMarkers_Full_",sampleName,".xlsx"))

########################################################################################################################

#split.by argument splits featureplot by sample!! (1-2-3)
FeaturePlot(object = seuratObjNew, features = "Apoe", cols = c("grey", "blue"),reduction = "SCT_umap", 
            min.cutoff = 'q2', max.cutoff = 'q98', split.by = "hash.ID", pt.size = 1.5, order = T) 
#split.by argument splits featureplot by genotype!! (DKO - WT)
FeaturePlot(object = seuratObjNew, features = "Apoe", cols = c("grey", "blue"),reduction = "SCT_umap", 
            min.cutoff = 'q2', max.cutoff = 'q98', split.by = "Genotype", pt.size = 1.5, order = T) 

FeaturePlot(object = seuratObjNew, features = "Apol7c", cols = c("grey", "blue"),reduction = "SCT_umap", 
            min.cutoff = 'q2', max.cutoff = 'q98', split.by = "Genotype", pt.size = 1.5, order = T) 

FeaturePlot(object = seuratObjNew, features = "Vim", cols = c("grey", "blue"),reduction = "SCT_umap", 
            min.cutoff = 'q2', max.cutoff = 'q98', split.by = "Genotype", pt.size = 1.5, order = T) 

####################################################################################################################

## Check important markers!!
F1<-FeaturePlot(object = seuratObjNew, features = c("Itgae","CD103","CD207-mh","CD62L"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5,order = T)
ggsave(F1, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/4_Featureplot_for_new_annotation",sampleName,"_2021.png"), height = 14, width = 20, dpi = "retina")

## Create final annotation of subset !!Happy with minim annot. -> Keep annotation!!
seuratObjNew@meta.data$annotated_clusters_final2021 <- seuratObjNew@meta.data$annotated_clusters_minim

U4 <- DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters_final2021", label.size = 4, cols = Colorset)
ggsave(U4, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/5_UMAP_annotated_subset_",sampleName,"_final2021.png"), height = 10, width = 15, dpi = "retina")

U5<-DimPlot(seuratObjNew, reduction = "SCT_umap", label = F, group.by = "Genotype", label.size = 5, cols = c("Green","Orange"))
ggsave(U5, file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/5_UMAP_spleen_local_",sampleName,"_final2021.png"), height = 10, width = 15, dpi = "retina")

########################################################################################################################

Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters_final2021

# create a dataset
Sample <- seuratObjNew@meta.data$orig.ident
cluster <- seuratObjNew@active.ident
Aggr <- rep(sampleName,length(cluster)) 

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))

# Stacked
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")


##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="SAM05"),]
data_split2<-data[which(data$Sample=="SAM06"),]
data_split3<-data[which(data$Sample=="SAM07"),]
data_split4<-data[which(data$Sample=="SAM08"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4)

png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/6_SampleDistribution_ggplot2_cDC1s_",sampleName,"_adjusted_2021.png"), width = 2000, height = 1500, res = 300)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

############################################################################################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObjNew, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_cDC1s_",sampleName,"_final2021.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterfinal2021']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_cDC1s_clean_",sampleName,"_2021.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){

  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_SCTclus_cDC1s_",sampleName,"_final2021.xlsx"))

######################################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_SCTclus <- FindAllMarkers(seuratObjNew, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_SCTclus_cDC1s",sampleName,"_final2021.rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTclusterfinal2021']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_cDC1s_clean_",sampleName,"_2021.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-names(table(ADTMarkers_SCTclus$cluster))
ADTmarkersList_SCTclus<-list()

for(i in totalNrADTclusters_SCTclus){

  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_SCTclus_cDC1s_",sampleName,"_final2021.xlsx"))


########################################################################################################################

##### Read final object
seuratObjNew <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_cDC1s_clean_",sampleName,"_2021.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_cDC1s_clean_",sampleName,"_2021.rds"))

##### Save object
saveRDS(seuratObjNew, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_cDC1s_clean_",sampleName,"_2021.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_cDC1s_clean_",sampleName,"_2021.rds"))

########################################################################################################################
########################################################################################################################
########################################################################################################################

## Analysis paper 2023
## Use HTO information for paper Victor (17/02/23)
## Need to filter the object: remove doublet and negative MULTI_ID cells!!! New object!!

table(seuratObjNew@meta.data$MULTI_ID)
# Doublet Hashtag1 Hashtag2 Hashtag3 Hashtag4 Hashtag5 Hashtag6 Negative 
# 1527     4043     2738     3868     4102     3748     4527       94

Idents(seuratObjNew)<-seuratObjNew@meta.data$MULTI_ID
levels(Idents(seuratObjNew))

seuratObjNew2023<-subset(seuratObjNew, idents = levels(Idents(seuratObjNew))[-c(4,5)])

seuratObjNew2023@meta.data$replicate_info<-paste0(seuratObjNew2023@meta.data$orig.ident,"_",seuratObjNew2023@meta.data$MULTI_ID)
levels(as.factor(seuratObjNew2023@meta.data$replicate_info))

seuratObjNew2023@meta.data$annotated_clusters_final2021_v4<-seuratObjNew2023@meta.data$annotated_clusters_final2021
levels(seuratObjNew2023@meta.data$annotated_clusters_final2021_v4)[1]<-"Pre-cDC1s"
levels(seuratObjNew2023@meta.data$annotated_clusters_final2021_v4)[5]<-"Cxcl9+ cDC1s"
Idents(seuratObjNew2023)<-seuratObjNew2023@meta.data$annotated_clusters_final2021_v4

#Save as pdf
Colorset<-c(brewer.pal(12,"Set3"))

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/7_Figure_cDC1s_clean_",sampleName,"_paper_2023.pdf"), height = 10, width = 15)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters_final2021_v4", label.size = 5,
        cols =c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) 
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/7_Figure_cDC1s_clean_",sampleName,"_split_paper_2023.pdf"), height = 10, width = 27)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", label = F, repel = T, group.by = "annotated_clusters_final2021_v4", label.size = 5, pt.size = 0.7,
        split.by = "Genotype", 
        cols =c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) 
dev.off()

######################################################################

## Check clustering for Victor: Early <-> Late immature cDC1s
## Perhaps Harmony not working optimally??

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator

FeatureScatter(seuratObjNew2023, feature1 = "adt_CD62L", feature2 = "adt_CD43", cols = Colorset, group.by = "annotated_clusters_final2021_v4")

Idents(seuratObjNew2023)<-seuratObjNew2023$Genotype
seuratObjNew2023_DKO<-subset(seuratObjNew2023, ident = "DKO")
seuratObjNew2023_WT<-subset(seuratObjNew2023, ident = "WT")

Idents(seuratObjNew2023_DKO)<-seuratObjNew2023_DKO$annotated_clusters_final2021_v4
Idents(seuratObjNew2023_WT)<-seuratObjNew2023_WT$annotated_clusters_final2021_v4

for (i in 1:length(levels(seuratObjNew2023$annotated_clusters_final2021_v4))){
  Cells_DKO<-rownames(seuratObjNew2023_DKO@meta.data[which(seuratObjNew2023_DKO$annotated_clusters_final2021_v4 == levels(seuratObjNew2023_DKO$annotated_clusters_final2021_v4)[i]),])
  F1_DKO<-FeatureScatter(seuratObjNew2023_DKO, feature1 = "adt_CD62L", feature2 = "adt_CD43", cols = Colorset[i], cells = Cells_DKO) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_DKO"))
  F2_DKO<-FeatureScatter(seuratObjNew2023_DKO, feature1 = "adt_CD102", feature2 = "adt_CD22", cols = Colorset[i], cells = Cells_DKO) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_DKO"))
  F3_DKO<-FeatureScatter(seuratObjNew2023_DKO, feature1 = "adt_CD103", feature2 = "adt_CD207-mh", cols = Colorset[i], cells = Cells_DKO) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_DKO"))
  
  Cells_WT<-rownames(seuratObjNew2023_WT@meta.data[which(seuratObjNew2023_WT$annotated_clusters_final2021_v4 == levels(seuratObjNew2023_WT$annotated_clusters_final2021_v4)[i]),])
  F1_WT<-FeatureScatter(seuratObjNew2023_WT, feature1 = "adt_CD62L", feature2 = "adt_CD43", cols = Colorset[i], cells = Cells_WT) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_WT"))
  F2_WT<-FeatureScatter(seuratObjNew2023_WT, feature1 = "adt_CD102", feature2 = "adt_CD22", cols = Colorset[i], cells = Cells_WT) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_WT"))
  F3_WT<-FeatureScatter(seuratObjNew2023_WT, feature1 = "adt_CD103", feature2 = "adt_CD207-mh", cols = Colorset[i], cells = Cells_WT) +
    xlim(0,5) + ylim(0,5) + ggtitle(paste0(levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_WT"))
  
  pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Scatterplots_for_Victor_20022023_",levels(seuratObjNew2023$annotated_clusters_final2021_v4)[i],"_",sampleName,".pdf"), height = 7, width = 12)
  print(F1_DKO+F1_WT)
  print(F2_DKO+F2_WT)
  print(F3_DKO+F3_WT)
  dev.off()
}

## Check congruency between annotations
## Read harmony object
seuratObjNew2023_Harmony <- readRDS(file=paste0(sampleFolder,"results_harmony/Robjects/seuratObj_paper_SAM_merge_harmony_2023.rds"))
intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)) #22993
setdiff(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)) #33
setdiff(colnames(seuratObjNew2023_Harmony),colnames(seuratObjNew2023)) #0

##Subset to same cells
seuratObjNew2023_Harmony<-subset(seuratObjNew2023_Harmony, cells = intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)))
seuratObjNew2023<-subset(seuratObjNew2023, cells = intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)))

## Check frequency tables
table(seuratObjNew2023$annotated_clusters_final2021_v4)
table(seuratObjNew2023_Harmony$annotated_clusters_final2021_v4)

## Check sankey plot between them
Merge_meta<-seuratObjNew2023@meta.data
Harmony_meta<-seuratObjNew2023_Harmony@meta.data

all(rownames(Merge_meta) == rownames(Harmony_meta))
Test<-as.data.frame(cbind(as.character(Merge_meta$annotated_clusters_final2021_v4),
                          as.character(Harmony_meta$annotated_clusters_final2021_v4)))
rownames(Test)<-rownames(Merge_meta)
colnames(Test)<-c("Annotation_no_Harmony","Annotation_Harmony")

# select your annotations you want to correspond
annots <- dplyr::select(Test,Annotation_no_Harmony,Annotation_Harmony) 

# summarise them with count
annot.tab <- annots %>%
  dplyr::count(Annotation_no_Harmony,Annotation_Harmony) 

colnames(annot.tab) <- c("Non-harmony", "Harmony", "value")
annot.tab$Harmony <- paste(annot.tab$Harmony, " ", sep="")

# # From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(annot.tab$`Non-harmony`), annot.tab$Harmony)) %>% unique()

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
annot.tab$IDnonHarmony=match(annot.tab$`Non-harmony`, nodes$name)-1 
annot.tab$IDHarmony=match(annot.tab$Harmony, nodes$name)-1

# Add a 'group' column to each connection:
annot.tab$group <- as.factor(c(rep("type_a",7),rep("type_b",1),rep("type_a",29)))

# Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c("Cxcl9","Early_imm","Early_mat","Late_imm","Late_mat","Pre","Prolif",
                           "Cxcl9","Early_imm","Late_imm","Prolif","Early_mat","Pre","Late_mat"))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "Cxcl9","Early_imm","Early_mat","Late_imm","Late_mat","Pre","Prolif"]) .range(["gray", "gold", "#80B1D3","#BEBADA","#FDB462","#FB8072","#B3DE69","#8DD3C7","#FFFFB3"])'

# Make the Network
p <- sankeyNetwork(Links = annot.tab, Nodes = nodes,
                   Source = "IDnonHarmony", Target = "IDHarmony",
                   Value = "value", NodeID = "name", 
                   colourScale=my_color, LinkGroup="group", NodeGroup="group",
                   sinksRight=FALSE, nodeWidth=30, fontSize=10, nodePadding=15) #40, 15, 20 iterations = 0

p

# save the widget
saveWidget(p, file=paste0( sampleFolder, "results_merge_non_harmony/Annotation/8_Comparison_Harmony_nonHarmony_sankeyColor.html"))

############################################################################################################

## Redo annotation: check higher resolution clustering
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1", repel = T, label = T, label.size = 2)

FeaturePlot(object = seuratObjNew2023, features = c("Nr4a3","Cd207","Itgae"), cols = c("grey", "blue"),reduction = "SCT_umap", 
            min.cutoff = 'q2', max.cutoff = 'q98', split.by = "Genotype", pt.size = 1.5, order = T, ncol = 3) 

## Redo clustering with 1.2/1.5 res!!
Perplexity<-30

seuratObjNew2023 <- FindNeighbors(object = seuratObjNew2023, reduction = "SCT_pca", dims = 1:Perplexity)
resolutions<-c(1.2,1.3,1.4,1.5)
for(res in resolutions){
  seuratObjNew2023 <- FindClusters(object = seuratObjNew2023,  graph.name = "SCT_snn", resolution = res)
}

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9.Overview_clustering_UMAP_",sampleName,".pdf"), height = 7, width = 10)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "annotated_clusters_final2021_v4", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "annotated_clusters", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.0.8", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1.2", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1.3", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1.4", repel = T, label = T, label.size = 4)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1.5", repel = T, label = T, label.size = 4)
dev.off()

## Decide to continue with 1.4
DimPlot(seuratObjNew2023, reduction = "SCT_umap", group.by = "SCT_snn_res.1.4", repel = T, label = T, label.size = 4)
seuratObjNew2023$annotation_paper_2023_detailed<-seuratObjNew2023$SCT_snn_res.1.4
seuratObjNew2023$annotation_paper_2023<-seuratObjNew2023$SCT_snn_res.1.4
Idents(seuratObjNew2023)<-seuratObjNew2023$SCT_snn_res.1.4

## Check clusters
clusterMatrix<-seuratObjNew2023@meta.data
tsneTable<-as.data.frame(seuratObjNew2023[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObjNew2023[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_annotation_Color_SCT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObjNew2023@meta.data$SCT_snn_res.1.4))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObjNew2023, cells = rownames(seuratObjNew2023@meta.data[which(seuratObjNew2023@meta.data$SCT_snn_res.1.4==i),])))
  C1<-C1+ggtitle(paste0("SCT_cluster_",i))
  print(C1)
}
dev.off()

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObjNew2023)<-seuratObjNew2023@meta.data$SCT_snn_res.1.4

ADTMarkers_SCTclus <- FindAllMarkers(seuratObjNew2023, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_SCTclus_numbered_paper_2023_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster2023']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-max(as.numeric(names(table(ADTMarkers_SCTclus$cluster))))
totalNrADTclusters_SCTclusPlusOne<-totalNrADTclusters_SCTclus+1
ADTmarkersList_SCTclus<-list()

for(i in 1:totalNrADTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrADTclusters_SCTclus)

### Write to Excel
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_SCTclus_numbered_paper_2023_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObjNew2023, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_numbered_paper_2023_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster2023']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-max(as.numeric(names(table(SCTMarkers_SCTclus$cluster))))
totalNrSCTclusters_SCTclusPlusOne<-totalNrSCTclusters_SCTclus+1
SCTmarkersList_SCTclus<-list()

for(i in 1:totalNrSCTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(SCTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrSCTclusters_SCTclus)

### Write to Excel
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_SCTclus_numbered_paper_2023_",sampleName,".xlsx"))

## New annotation
levels(seuratObjNew2023@meta.data$annotation_paper_2023_detailed) <- c("Early immature cDC1s 1","Nr4a3+CD207+ Late imm cDC1s","Egr+Nr4a+ Late imm cDC1s 1",
                                                                "Egr+Nr4a+ Late imm cDC1s 2","Prolif cDC1s 4","Prolif cDC1s 3",
                                                                "Apoe+Itgae+CD207+ Late imm cDC1s","CD62L+ Early imm cDC1s 2","Late mat cDC1s",
                                                                "Hspa+ Early imm cDC1s 3","Apoe+ Late imm cDC1s 1","Apoe+ Late imm cDC1s 2",
                                                                "Prolif cDC1s 1","Early mat cDC1s","Pre-cDC1s","Cxcl9+ cDC1s 1",
                                                                "Hspa+ Early imm cDC1s 4","Cxcl9+ cDC1s 2","Prolif cDC1s 2",
                                                                "CD62L+Cxcl9+ cDC1s 3")
levels(seuratObjNew2023@meta.data$annotation_paper_2023) <- c("Early immature cDC1s","Late immature cDC1s","Late immature cDC1s",
                                                              "Late immature cDC1s","Proliferating cDC1s","Proliferating cDC1s",
                                                              "Late immature cDC1s","Early immature cDC1s","Late mature cDC1s",
                                                              "Early immature cDC1s","Late immature cDC1s","Late immature cDC1s",
                                                              "Proliferating cDC1s","Early mature cDC1s","Pre-cDC1s","Cxcl9+ cDC1s",
                                                              "Early immature cDC1s","Cxcl9+ cDC1s","Proliferating cDC1s",
                                                              "Cxcl9+ cDC1s")
seuratObjNew2023@meta.data$annotation_paper_2023<-factor(seuratObjNew2023@meta.data$annotation_paper_2023,
                                                        as.character(c("Pre-cDC1s","Proliferating cDC1s","Early immature cDC1s","Late immature cDC1s",
                                                                       "Cxcl9+ cDC1s","Early mature cDC1s",'Late mature cDC1s'))) #reorder levels

Idents(seuratObjNew2023)<-seuratObjNew2023@meta.data$annotation_paper_2023

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_Figure_new_annotation_detailed_",sampleName,"_paper_2023.pdf"), height = 10, width = 15)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", label = T, repel = T, group.by = "annotation_paper_2023_detailed", label.size = 5) 
dev.off()

Colorset<-c(brewer.pal(12,"Set3"))

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_Figure_new_annotation_",sampleName,"_paper_2023.pdf"), height = 10, width = 15)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", label = T, repel = T, group.by = "annotation_paper_2023", label.size = 5,
        cols =c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) 
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_Figure_new_annotation_",sampleName,"_split_paper_2023.pdf"), height = 10, width = 27)
DimPlot(seuratObjNew2023, reduction = "SCT_umap", label = F, repel = T, group.by = "annotation_paper_2023", label.size = 5, pt.size = 0.7,
        split.by = "Genotype", 
        cols =c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) 
dev.off()

################################################
################################################

# create a dataset
Sample <- seuratObjNew2023@meta.data$replicate_info
cluster <- seuratObjNew2023@active.ident
Aggr <- rep(sampleName,length(cluster)) 

data <- data.frame(table(Sample, cluster))
data2 <- data.frame(table(cluster,Aggr))
data4<-data.frame(table(cluster,Sample)) #Added 27/07/23

# Stacked
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")

png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_SampleDistribution_ggplot2_cDC1s_",sampleName,"_paper_2023.png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('cadetblue1',"cadetblue2",'cadetblue3',
                                                                                            'royalblue2',"royalblue3",'royalblue4',
                                                                                            'indianred1','indianred2','indianred3',
                                                                                            'firebrick2','firebrick3','firebrick4')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
## Removed round function!!! (03/03/23)
data_split1<-data[which(data$Sample=="SAM05_Hashtag1"),]
data_split2<-data[which(data$Sample=="SAM05_Hashtag2"),]
data_split3<-data[which(data$Sample=="SAM05_Hashtag3"),]
data_split4<-data[which(data$Sample=="SAM06_Hashtag1"),]
data_split5<-data[which(data$Sample=="SAM06_Hashtag2"),]
data_split6<-data[which(data$Sample=="SAM06_Hashtag3"),]
data_split7<-data[which(data$Sample=="SAM07_Hashtag4"),]
data_split8<-data[which(data$Sample=="SAM07_Hashtag5"),]
data_split9<-data[which(data$Sample=="SAM07_Hashtag6"),]
data_split10<-data[which(data$Sample=="SAM08_Hashtag4"),]
data_split11<-data[which(data$Sample=="SAM08_Hashtag5"),]
data_split12<-data[which(data$Sample=="SAM08_Hashtag6"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-(x/sum(data_split1$Freq))*100})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-(x/sum(data_split2$Freq))*100})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-(x/sum(data_split3$Freq))*100})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-(x/sum(data_split4$Freq))*100})
data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-(x/sum(data_split5$Freq))*100})
data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-(x/sum(data_split6$Freq))*100})
data_split7$Freq<-lapply(data_split7$Freq,function(x){x<-(x/sum(data_split7$Freq))*100})
data_split8$Freq<-lapply(data_split8$Freq,function(x){x<-(x/sum(data_split8$Freq))*100})
data_split9$Freq<-lapply(data_split9$Freq,function(x){x<-(x/sum(data_split9$Freq))*100})
data_split10$Freq<-lapply(data_split10$Freq,function(x){x<-(x/sum(data_split10$Freq))*100})
data_split11$Freq<-lapply(data_split11$Freq,function(x){x<-(x/sum(data_split11$Freq))*100})
data_split12$Freq<-lapply(data_split12$Freq,function(x){x<-(x/sum(data_split12$Freq))*100})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4,data_split5,data_split6,data_split7,data_split8,
                data_split9,data_split10,data_split11,data_split12)

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_SampleDistribution_ggplot2_cDC1s_",sampleName,"_adjusted_paper_2023_correct.pdf"), width = 12, height = 9)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('cadetblue1',"cadetblue2",'cadetblue3',
                                                                                                'royalblue2',"royalblue3",'royalblue4',
                                                                                                'indianred1','indianred2','indianred3',
                                                                                                'firebrick2','firebrick3','firebrick4')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

## Extra variation -> combine after correcting like this (RP and WP separate!!)
data_new$Genotype<-data_new$Sample
data_new$Genotype<-gsub("SAM05","WT",data_new$Genotype)
data_new$Genotype<-gsub("SAM06","WT",data_new$Genotype)
data_new$Genotype<-gsub("SAM07","DKO",data_new$Genotype)
data_new$Genotype<-gsub("SAM08","DKO",data_new$Genotype)

data_new$Freq<-as.numeric(data_new$Freq)
data_new_combo<-aggregate(Freq ~ Genotype + cluster, data = data_new, FUN = sum, na.rm = TRUE)
data_new_combo$Genotype<-factor(data_new_combo$Genotype, levels = c("WT_Hashtag1", "WT_Hashtag2", "WT_Hashtag3", 
                                                                    "DKO_Hashtag4", "DKO_Hashtag5", "DKO_Hashtag6"))

## Divide all values by 2 to create value out of 100 again per celltype/replicate. True percentages...
data_new_combo$Freq<-data_new_combo$Freq/2

## Try example with 10k cells -> *100 and round
data_new_combo_10k<-data_new_combo
data_new_combo_10k$Freq<-round(data_new_combo_10k$Freq*100,0)

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_SampleDistribution_ggplot2_cDC1s_v3_",sampleName,"_adjusted_paper_2023_correct.pdf"), width = 12, height = 9)
ggplot(data_new_combo, aes(fill=Genotype, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('royalblue2',"royalblue3",'royalblue4',
                                                                                                        'firebrick2','firebrick3','firebrick4')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

## Convert and save data tables
data_wide<-reshape(data, idvar = "Sample", timevar = "cluster", direction = "wide")
data_new_wide<-reshape(data_new[-4], idvar = "Sample", timevar = "cluster", direction = "wide")
data_new_combo_wide<-reshape(data_new_combo, idvar = "Genotype", timevar = "cluster", direction = "wide")
data_new_combo_10k_wide<-reshape(data_new_combo_10k, idvar = "Genotype", timevar = "cluster", direction = "wide")

write.xlsx(data_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table1_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table2_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_combo_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table3_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_combo_10k_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table3_10k_",sampleName,"_paper_2023.xlsx"))

## Check with independent t-test (https://www.scribbr.com/statistics/statistical-tests/)
## Check per cluster
data_new_combo1<-data_new_combo[which(data_new_combo$cluster == "Pre-cDC1s"),c(1,3)]
data_new_combo1$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo1$Genotype)
data_new_combo2<-data_new_combo[which(data_new_combo$cluster == "Proliferating cDC1s"),c(1,3)]
data_new_combo2$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo2$Genotype)
data_new_combo3<-data_new_combo[which(data_new_combo$cluster == "Early immature cDC1s"),c(1,3)]
data_new_combo3$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo3$Genotype)
data_new_combo4<-data_new_combo[which(data_new_combo$cluster == "Late immature cDC1s"),c(1,3)]
data_new_combo4$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo4$Genotype)
data_new_combo5<-data_new_combo[which(data_new_combo$cluster == "Cxcl9+ cDC1s"),c(1,3)]
data_new_combo5$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo5$Genotype)
data_new_combo6<-data_new_combo[which(data_new_combo$cluster == "Early mature cDC1s"),c(1,3)]
data_new_combo6$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo6$Genotype)
data_new_combo7<-data_new_combo[which(data_new_combo$cluster == "Late mature cDC1s"),c(1,3)]
data_new_combo7$Genotype<-gsub(paste0("_Hashtag","[0-9]"),"",data_new_combo7$Genotype)

res1<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo1)
res2<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo2)
res3<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo3)
res4<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo4)
res5<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo5)
res6<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo6)
res7<-t.test(Freq ~ Genotype, var.equal=TRUE, data = data_new_combo7)

Test<-cbind(rbind(res1$estimate,res2$estimate,res3$estimate,res4$estimate,res5$estimate,res6$estimate,res7$estimate),
      rbind(res1$p.value,res2$p.value,res3$p.value,res4$p.value,res5$p.value,res6$p.value,res7$p.value))
colnames(Test)[3]<-"p.value"
rownames(Test)<-levels(Idents(seuratObjNew2023))

write.xlsx(as.data.frame(Test), paste0(sampleFolder,"results_merge_non_harmony/Stat_test/T_test_results_table3_",sampleName,"_paper_2023.xlsx"),
           rowNames = T)

# ###################################
# ## Other version for Victor: Combine WP and RP for each genotype
# seuratObjNew2023@meta.data$replicate_info_v2<-paste0(seuratObjNew2023@meta.data$Genotype,"_",seuratObjNew2023@meta.data$MULTI_ID)
# seuratObjNew2023@meta.data$replicate_info_v2<-factor(seuratObjNew2023@meta.data$replicate_info_v2, 
#                                                      levels = c("WT_Hashtag1", "WT_Hashtag2", "WT_Hashtag3", 
#                                                                 "DKO_Hashtag4", "DKO_Hashtag5", "DKO_Hashtag6"))
# 
# # create a dataset
# Sample <- seuratObjNew2023@meta.data$replicate_info_v2
# cluster <- seuratObjNew2023@active.ident
# Aggr <- rep(sampleName,length(cluster)) 
# 
# data <- data.frame(table(Sample, cluster))
# data2 <- data.frame(table(cluster,Aggr))
# 
# # barplotAggr(seuratObjNew, listLabels)
# 
# # Stacked
# ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
#   geom_bar(position="fill", stat="identity", colour="white")
# 
# png(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/12_SampleDistribution_ggplot2_cDC1s_v2_",sampleName,"_paper_2023.png"), width = 2000, height = 1500, res = 300)
# ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('royalblue2',"royalblue3",'royalblue4',
#                                                                                             'firebrick2','firebrick3','firebrick4')) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
#   geom_bar(position="fill", stat="identity", colour="white")
# dev.off()
# 
# ##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
# data_split1<-data[which(data$Sample=="WT_Hashtag1"),]
# data_split2<-data[which(data$Sample=="WT_Hashtag2"),]
# data_split3<-data[which(data$Sample=="WT_Hashtag3"),]
# data_split4<-data[which(data$Sample=="DKO_Hashtag4"),]
# data_split5<-data[which(data$Sample=="DKO_Hashtag5"),]
# data_split6<-data[which(data$Sample=="DKO_Hashtag6"),]
# data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
# data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
# data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
# data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})
# data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-round((x/sum(data_split5$Freq))*100,2)})
# data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-round((x/sum(data_split6$Freq))*100,2)})
# 
# data_new<-rbind(data_split1,data_split2,data_split3,data_split4,data_split5,data_split6)
# 
# pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/12_SampleDistribution_ggplot2_cDC1s_v2_",sampleName,"_adjusted_paper_2023.pdf"), width = 12, height = 9)
# ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=c('royalblue2',"royalblue3",'royalblue4',
#                                                                                                 'firebrick2','firebrick3','firebrick4')) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
#   geom_bar(position="fill", stat="identity", colour="white")
# dev.off()
# 
# ###################################

## Extra test for Victor (added 27/07/23)
ggplot(data4, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white") + theme(legend.position= "none")

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
## Removed round function!!! (03/03/23)
data_split1<-data4[which(data4$Sample=="SAM05_Hashtag1"),]
data_split2<-data4[which(data4$Sample=="SAM05_Hashtag2"),]
data_split3<-data4[which(data4$Sample=="SAM05_Hashtag3"),]
data_split4<-data4[which(data4$Sample=="SAM06_Hashtag1"),]
data_split5<-data4[which(data4$Sample=="SAM06_Hashtag2"),]
data_split6<-data4[which(data4$Sample=="SAM06_Hashtag3"),]
data_split7<-data4[which(data4$Sample=="SAM07_Hashtag4"),]
data_split8<-data4[which(data4$Sample=="SAM07_Hashtag5"),]
data_split9<-data4[which(data4$Sample=="SAM07_Hashtag6"),]
data_split10<-data4[which(data4$Sample=="SAM08_Hashtag4"),]
data_split11<-data4[which(data4$Sample=="SAM08_Hashtag5"),]
data_split12<-data4[which(data4$Sample=="SAM08_Hashtag6"),]

data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-(x/sum(data_split1$Freq))*100})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-(x/sum(data_split2$Freq))*100})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-(x/sum(data_split3$Freq))*100})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-(x/sum(data_split4$Freq))*100})
data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-(x/sum(data_split5$Freq))*100})
data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-(x/sum(data_split6$Freq))*100})
data_split7$Freq<-lapply(data_split7$Freq,function(x){x<-(x/sum(data_split7$Freq))*100})
data_split8$Freq<-lapply(data_split8$Freq,function(x){x<-(x/sum(data_split8$Freq))*100})
data_split9$Freq<-lapply(data_split9$Freq,function(x){x<-(x/sum(data_split9$Freq))*100})
data_split10$Freq<-lapply(data_split10$Freq,function(x){x<-(x/sum(data_split10$Freq))*100})
data_split11$Freq<-lapply(data_split11$Freq,function(x){x<-(x/sum(data_split11$Freq))*100})
data_split12$Freq<-lapply(data_split12$Freq,function(x){x<-(x/sum(data_split12$Freq))*100})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4,data_split5,data_split6,data_split7,data_split8,
                data_split9,data_split10,data_split11,data_split12)

ggplot(data_new, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() + 
  scale_fill_manual(values=c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")

## Extra variation -> combine after correcting like this (RP and WP separate!!)
data_new$Genotype<-data_new$Sample
data_new$Genotype<-gsub("SAM05","WT",data_new$Genotype)
data_new$Genotype<-gsub("SAM06","WT",data_new$Genotype)
data_new$Genotype<-gsub("SAM07","DKO",data_new$Genotype)
data_new$Genotype<-gsub("SAM08","DKO",data_new$Genotype)
data_new$Genotype2<-data_new$Sample
data_new$Genotype2<-gsub("SAM05_Hashtag1","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM05_Hashtag2","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM05_Hashtag3","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM06_Hashtag1","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM06_Hashtag2","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM06_Hashtag3","WT",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM07_Hashtag4","DKO",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM07_Hashtag5","DKO",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM07_Hashtag6","DKO",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM08_Hashtag4","DKO",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM08_Hashtag5","DKO",data_new$Genotype2)
data_new$Genotype2<-gsub("SAM08_Hashtag6","DKO",data_new$Genotype2)
data_new$cluster2<-paste0(data_new$Genotype,"_",data_new$cluster)

data_new$Freq<-as.numeric(data_new$Freq)
data_new_combo<-aggregate(Freq ~ cluster2 + Genotype2, data = data_new, FUN = sum, na.rm = TRUE)
data_new_combo$Genotype2<-factor(data_new_combo$Genotype2, levels = c("WT", "DKO"))
data_new_combo$cluster2<-factor(data_new_combo$cluster2, levels = levels(as.factor(data_new_combo$cluster2))[c(27,34,41,
                                                                                                              28,35,42,
                                                                                                              23,30,37,
                                                                                                              25,32,39,
                                                                                                              22,29,36,
                                                                                                              24,31,38,
                                                                                                              26,33,40,
                                                                                                              6,13,20,
                                                                                                              7,14,21,
                                                                                                              2,9,16,
                                                                                                              4,11,18,
                                                                                                              1,8,15,
                                                                                                              3,10,17,
                                                                                                              5,12,19)])

## Divide all values by 6 to create value out of 100 again per celltype/replicate. True percentages...
data_new_combo$Freq<-data_new_combo$Freq/6

## Try example with 10k cells -> *100 and round
data_new_combo_10k<-data_new_combo
data_new_combo_10k$Freq<-round(data_new_combo_10k$Freq*100,0)

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_SampleDistribution_ggplot2_cDC1s_v4_",sampleName,"_adjusted_paper_2023_correct.pdf"), width = 15, height = 9)
ggplot(data_new_combo, aes(fill=cluster2, y=Freq, x=Genotype2)) + theme_bw() + 
  scale_fill_manual(values=rep(c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)]), times = 2, each = 3)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

## Convert and save data tables
data_wide<-reshape(data, idvar = "cluster", timevar = "Sample", direction = "wide")
data_new_wide<-reshape(data_new[-c(4,5,6)], idvar = "cluster", timevar = "Sample", direction = "wide")
data_new_combo_wide<-reshape(data_new_combo, idvar = "cluster2", timevar = "Genotype2", direction = "wide")
data_new_combo_10k_wide<-reshape(data_new_combo_10k, idvar = "cluster2", timevar = "Genotype2", direction = "wide")

write.xlsx(data_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table1_v4_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table2_v4_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_combo_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table3_v4_",sampleName,"_paper_2023.xlsx"))
write.xlsx(data_new_combo_10k_wide, paste0(sampleFolder,"results_merge_non_harmony/Stat_test/Frequency_table3_10k_v4_",sampleName,"_paper_2023.xlsx"))

########################################################################################################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObjNew2023, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/SCTmarkersList_SCTclus_cDC1s_",sampleName,"_paper_2023.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterpaper_2023']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){

  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/SCTmarkersList_SCTclus_cDC1s_",sampleName,"_paper_2023.xlsx"))

######################################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_SCTclus <- FindAllMarkers(seuratObjNew2023, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/ADTmarkersList_SCTclus_cDC1s_",sampleName,"_paper_2023.rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTclusterpaper_2023']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-names(table(ADTMarkers_SCTclus$cluster))
ADTmarkersList_SCTclus<-list()

for(i in totalNrADTclusters_SCTclus){

  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_merge_non_harmony/Marker_lists/ADTmarkersList_SCTclus_cDC1s_",sampleName,"_paper_2023.xlsx"))

##################################################################################################################################

#################################################################
########## GET DIFF MARKERS DKO vs WT ##########
#################################################################

### Create new clusters: split on source
seuratObjNew2023@meta.data$newClustersTmp<-Idents(seuratObjNew2023)
seuratObjNew2023@meta.data$newClusters_paper<-paste0(seuratObjNew2023@meta.data$newClustersTmp,"_",seuratObjNew2023@meta.data$Genotype)
head(seuratObjNew2023@meta.data)

### Use the new clusters
seuratObjNew2023@meta.data$newClusters_paper<- as.factor(seuratObjNew2023@meta.data$newClusters_paper) #reorder levels
Idents(seuratObjNew2023)<-seuratObjNew2023@meta.data$newClusters_paper

########## 2. GET MARKERS  ##########
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew2023, ident.1 = ident1, ident.2 = ident2,
                             min.pct = 0.10) # No min diff pct!! , min.diff.pct = 0.15 or logFC logfc.threshold = 0.30,
  markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]

  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001

  markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
                     markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}

# ########## 2. GET MARKERS (everything!!) ##########
# getDEgenes<-function(ident1, ident2){
#   markersDiff <- FindMarkers(seuratObjNew2023, ident.1 = ident1, ident.2 = ident2,
#                              logfc.threshold = 0, min.pct = 0.01) #0.30, 0.10, min.diff.pct = 0.15
#   #markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
#   markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]
#   
#   markersDiff$geneSymbol<-rownames(markersDiff)
#   markersDiff$pct.1<-markersDiff$pct.1+0.001
#   markersDiff$pct.2<-markersDiff$pct.2+0.001
#   
#   markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
#                      markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
#   markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
#   return(markersDiff)
# }

#### Get diff markers WP vs RP #####
levels(Idents(seuratObjNew2023))

PrecDC1s_DKOvsWT<-getDEgenes("Pre-cDC1s_DKO","Pre-cDC1s_WT")
PrecDC1s_DKOvsWT<-PrecDC1s_DKOvsWT[order(PrecDC1s_DKOvsWT$avg_logFC,decreasing = T),]
head(PrecDC1s_DKOvsWT)
dim(PrecDC1s_DKOvsWT)

Proliferating_cDC1s_DKOvsWT<-getDEgenes("Proliferating cDC1s_DKO","Proliferating cDC1s_WT")
Proliferating_cDC1s_DKOvsWT<-Proliferating_cDC1s_DKOvsWT[order(Proliferating_cDC1s_DKOvsWT$avg_logFC,decreasing = T),]
head(Proliferating_cDC1s_DKOvsWT)
dim(Proliferating_cDC1s_DKOvsWT)

Early_imm_DKOvsWT<-getDEgenes("Early immature cDC1s_DKO","Early immature cDC1s_WT")
Early_imm_DKOvsWT<-Early_imm_DKOvsWT[order(Early_imm_DKOvsWT$avg_logFC,decreasing = T),]
head(Early_imm_DKOvsWT)
dim(Early_imm_DKOvsWT)

Late_imm_DKOvsWT<-getDEgenes("Late immature cDC1s_DKO","Late immature cDC1s_WT")
Late_imm_DKOvsWT<-Late_imm_DKOvsWT[order(Late_imm_DKOvsWT$avg_logFC,decreasing = T),]
head(Late_imm_DKOvsWT)
dim(Late_imm_DKOvsWT)

Cxcl9_cDC1s_DKOvsWT<-getDEgenes("Cxcl9+ cDC1s_DKO","Cxcl9+ cDC1s_WT")
Cxcl9_cDC1s_DKOvsWT<-Cxcl9_cDC1s_DKOvsWT[order(Cxcl9_cDC1s_DKOvsWT$avg_logFC,decreasing = T),]
head(Cxcl9_cDC1s_DKOvsWT)
dim(Cxcl9_cDC1s_DKOvsWT)

Early_mat_DKOvsWT<-getDEgenes("Early mature cDC1s_DKO","Early mature cDC1s_WT")
Early_mat_DKOvsWT<-Early_mat_DKOvsWT[order(Early_mat_DKOvsWT$avg_logFC,decreasing = T),]
head(Early_mat_DKOvsWT)
dim(Early_mat_DKOvsWT)

Late_mat_DKOvsWT<-getDEgenes("Late mature cDC1s_DKO","Late mature cDC1s_WT")
Late_mat_DKOvsWT<-Late_mat_DKOvsWT[order(Late_mat_DKOvsWT$avg_logFC,decreasing = T),]
head(Late_mat_DKOvsWT)
dim(Late_mat_DKOvsWT)

##add to list
listDiffMarkers<-tibble::lst(PrecDC1s_DKOvsWT,Proliferating_cDC1s_DKOvsWT,Early_imm_DKOvsWT,Late_imm_DKOvsWT,
                             Cxcl9_cDC1s_DKOvsWT,Early_mat_DKOvsWT,Late_mat_DKOvsWT)

lapply(listDiffMarkers, dim)
listDiffMarkers<-lapply(listDiffMarkers,function(x){x<-x[order(x$score, decreasing=T),]})

#Check settings
saveRDS(listDiffMarkers,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/markersDiffSamples_Full_",sampleName,"_paper_2023.rds"))

### Write to Excel
write.xlsx(listDiffMarkers, file = paste0(sampleFolder,"results_merge_non_harmony/Marker_lists/summaryDiffMarkers_Full_",sampleName,"_paper_2023.xlsx"))

########################################################################################################################

##### Read final object
seuratObjNew2023 <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_",sampleName,"_2023.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

##### Save paper object
saveRDS(seuratObjNew2023, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_",sampleName,"_2023.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

############################################################################################################

## Check congruency between annotations
## Read harmony object
seuratObjNew2023_Harmony <- readRDS(file=paste0(sampleFolder,"results_harmony/Robjects/seuratObj_paper_SAM_merge_harmony_2023.rds"))
intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)) #22993
setdiff(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)) #33
setdiff(colnames(seuratObjNew2023_Harmony),colnames(seuratObjNew2023)) #0

##Subset to same cells
seuratObjNew2023_Harmony<-subset(seuratObjNew2023_Harmony, cells = intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)))
seuratObjNew2023<-subset(seuratObjNew2023, cells = intersect(colnames(seuratObjNew2023),colnames(seuratObjNew2023_Harmony)))

## Check frequency tables
table(seuratObjNew2023$annotation_paper_2023)
table(seuratObjNew2023_Harmony$annotated_clusters_final2021_v4)

## Check sankey plot between them
Merge_meta<-seuratObjNew2023@meta.data
Harmony_meta<-seuratObjNew2023_Harmony@meta.data

all(rownames(Merge_meta) == rownames(Harmony_meta))
Test<-as.data.frame(cbind(as.character(Merge_meta$annotation_paper_2023),
                          as.character(Harmony_meta$annotated_clusters_final2021_v4)))
rownames(Test)<-rownames(Merge_meta)
colnames(Test)<-c("Annotation_no_Harmony","Annotation_Harmony")

# select your annotations you want to correspond
annots <- dplyr::select(Test,Annotation_no_Harmony,Annotation_Harmony) 

# summarise them with count
annot.tab <- annots %>%
  dplyr::count(Annotation_no_Harmony,Annotation_Harmony) 

colnames(annot.tab) <- c("Non-harmony", "Harmony", "value")
annot.tab$Harmony <- paste(annot.tab$Harmony, " ", sep="")

# # From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(annot.tab$`Non-harmony`), annot.tab$Harmony)) %>% unique()

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
annot.tab$IDnonHarmony=match(annot.tab$`Non-harmony`, nodes$name)-1 
annot.tab$IDHarmony=match(annot.tab$Harmony, nodes$name)-1

# Add a 'group' column to each connection:
annot.tab$group <- as.factor(c(rep("type_a",10),rep("type_b",1),rep("type_a",29)))

# Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c("Cxcl9","Early_imm","Early_mat","Late_imm","Late_mat","Pre","Prolif",
                           "Cxcl9","Early_imm","Early_mat","Late_imm","Late_mat","Pre","Prolif"))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "Cxcl9","Early_imm","Early_mat","Late_imm","Late_mat","Pre","Prolif"]) .range(["gray", "gold", "#80B1D3","#BEBADA","#FDB462","#FB8072","#B3DE69","#8DD3C7","#FFFFB3"])'

# Make the Network
p <- sankeyNetwork(Links = annot.tab, Nodes = nodes,
                   Source = "IDnonHarmony", Target = "IDHarmony",
                   Value = "value", NodeID = "name", 
                   colourScale=my_color, LinkGroup="group", NodeGroup="group",
                   sinksRight=FALSE, nodeWidth=30, fontSize=10, nodePadding=15) #40, 15, 20 iterations = 0

p


# save the widget
saveWidget(p, file=paste0( sampleFolder, "results_merge_non_harmony/Annotation/9_Comparison_Harmony_nonHarmony_sankeyColor_new_annotation.html"))

#################################################################################################

## Extra plots Victor (30/10/23)
seuratObj_filtered$Genotype<-as.factor(seuratObj_filtered$Genotype)
seuratObj_filtered$Genotype<-factor(seuratObj_filtered$Genotype, levels = c("WT","DKO"))

Colorset<-c(brewer.pal(12,"Set3"))

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Annotation/9_Figure_new_annotation_",sampleName,"_split_paper_30102023_different_order.pdf"), height = 10, width = 27)
DimPlot(seuratObj_filtered, reduction = "SCT_umap", label = F, repel = T, group.by = "annotation_paper_2023", label.size = 5, pt.size = 0.7,
        split.by = "Genotype", 
        cols =c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])) 
dev.off()

Extra_features<-c('Ccr7', 'Fscn1', 'Slco5a1', 'Tmem176a', 'Abcg1', 'Apoe', 'Apol7c', 'Apol10b',
                  'Atf4', 'Ddit3', 'Cars', 'Yars', 'Cox6a2', 'Dnase1l3', 'Il4i1', 'Mical3', 'H2-M2')

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Violinplots_for_Victor_30102023_split_by_genotype_different_order_v2_",sampleName,".pdf"), height = 7, width = 10)
for (feature in Extra_features) {
  V1<-VlnPlot(seuratObj_filtered, features = feature, split.by = "Genotype", cols = c("#ec672a","gray90"), pt.size = 0.3)
  print(V1)
}
dev.off()

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Feature_plots/Featureplots_for_Victor_30102023_split_by_genotype_different_order_",sampleName,"_blue_grey.pdf"), height = 7, width = 16)
for (feature in Extra_features) {
  F1<-FeaturePlot(object = seuratObj_filtered, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, 
                  split.by = "Genotype", order=T)
  print(F1)
}
dev.off()

################################################################################
########## FINAL OBJECT FOR PAPER AND ONLINE TOOL
################################################################################

## Filter ADTs for final object upload
## Remove ADTs which are superfluous! Suggestion from Niels!

## Load in ADT panel
Panel_Niels<-read.xlsx("CITE-seq_antibody_panel.xlsx")
Other_panel<-read.xlsx("PanelSAM_Clint.xlsx")
ADT_names<-rownames(seuratObjNew2023@assays$ADT)

intersect(Other_panel$AB_name,Panel_Niels$Target)
ADT_remove<-setdiff(Other_panel$AB_name,Panel_Niels$Target) #Missing 2! 19 instead of 21
setdiff(Panel_Niels$Target,Other_panel$AB_name)
setdiff(Panel_Niels$Target,intersect(Other_panel$AB_name,Panel_Niels$Target))

Test1<-sort(c(Panel_Niels$Target,rep("X",21)))
Test2<-sort(Other_panel$AB_name)

Testdf<-cbind(Test1,Test2)
Extra_ADT<-c("CD278.1","CD309-A0553") #2 missing ADT

rownames(Other_panel)<-make.unique(Other_panel$AB_name)
Other_panel_filtered<-Other_panel[!(rownames(Other_panel) %in% ADT_remove),] #Initial filter to 171

rownames(Other_panel_filtered)<-(Other_panel_filtered$whitelist_name)
Other_panel_filtered2<-Other_panel_filtered[!(rownames(Other_panel_filtered) %in% Extra_ADT),] #Subsequent filter to 169

intersect(Other_panel$AB_name,ADT_names)
intersect(Other_panel$whitelist_name,ADT_names)
intersect(Other_panel_filtered2$whitelist_name,ADT_names) #This is the list to filter on!!!! Matches seuratobj names!

seuratObj_filtered<-seuratObjNew2023
seuratObj_filtered@assays$ADT@counts<-seuratObj_filtered@assays$ADT@counts[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@data<-seuratObj_filtered@assays$ADT@data[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@scale.data<-seuratObj_filtered@assays$ADT@scale.data[intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name),]
seuratObj_filtered@assays$ADT@var.features<-intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name)
seuratObj_filtered@assays$ADT@meta.features<-seuratObj_filtered@assays$ADT@meta.features[Other_panel_filtered2$whitelist_name,]
# seuratObj_filtered@reductions$ADT_pca<-seuratObj_filtered@reductions$ADT_pca ##Can't subset pca!! Which dim to remove??

## Decided not to subset the metadata -> show the progression of the annotation
## The labels are clear enough for readers I think...

## Extra fix metadata
seuratObj_filtered$newClusters_paper<-factor(seuratObj_filtered$newClusters_paper, levels = levels(seuratObj_filtered$newClusters_paper)[c(11,12,13,14,3,4,7,8,1,2,5,6,9,10)])

########################################################################################################################

##### Read paper object
seuratObj_filtered <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_filtered_ADT_",sampleName,"_paper_2023.rds"))
diagnostics_filtered <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_filtered_ADT_",sampleName,"_paper_2023.rds"))

##### Save paper object
saveRDS(seuratObj_filtered, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_filtered_ADT_",sampleName,"_paper_2023.rds"))
saveRDS(diagnostics_filtered, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_filtered_ADT_",sampleName,"_paper_2023.rds"))

############################################################################################################

## Create diet object for online tool (01/08/23)
seuratObj_filtered_diet<-DietSeurat(seuratObj_filtered, counts = T, data = T, scale.data = F,
                                    assays = c("SCT","ADT"), dimreducs = "SCT_umap", graphs = NULL)

## New metadata names
## "orig.ident" "replicate_info" "Spleen_localisation" "Genotype" "SCT_snn_res.1.4" "annotation_paper_2023_detailed" "annotation_paper_2023" "newClusters_paper"
seuratObj_filtered_diet$Sample<-seuratObj_filtered_diet$orig.ident
seuratObj_filtered_diet$Replicate<-seuratObj_filtered_diet$replicate_info
seuratObj_filtered_diet$Numbered_SCT_clusters<-seuratObj_filtered_diet$SCT_snn_res.1.4
seuratObj_filtered_diet$Annotation_detailed<-seuratObj_filtered_diet$annotation_paper_2023_detailed
seuratObj_filtered_diet$Annotation<-seuratObj_filtered_diet$annotation_paper_2023
seuratObj_filtered_diet$Annotation_split_by_genotype<-seuratObj_filtered_diet$newClusters_paper

## Update detailed annotation Pre-cDC1s (same as annotation!!)
levels(seuratObj_filtered_diet$Annotation_detailed)[15]<-"Pre cDC1s"

DimPlot(seuratObj_filtered_diet, reduction = "SCT_umap", label = T, repel = T, group.by = "Annotation_detailed", label.size = 5,
        cols =gg_color_hue(20))

col2hex(c("#ec672a","gray90"))
c('cadetblue','royalblue','indianred','firebrick')
Colset<-c(brewer.pal(n = 7, name = "Set3")[1],"Yellow",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)])
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# ggplotColours <- function(n = 6, h = c(0, 360) + 15){ #Same results
#   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
All<-c(levels(as.factor(seuratObj_filtered_diet$Annotation)),
       levels(as.factor(seuratObj_filtered_diet$Sample)),levels(as.factor(seuratObj_filtered_diet$Replicate)),
       levels(as.factor(seuratObj_filtered_diet$Spleen_localisation)),levels(as.factor(seuratObj_filtered_diet$Genotype)),
       levels(as.factor(seuratObj_filtered_diet$Numbered_SCT_clusters)),levels(as.factor(seuratObj_filtered_diet$Annotation_detailed)),
       levels(as.factor(seuratObj_filtered_diet$Annotation_split_by_genotype)))
Color_info<-c(c(brewer.pal(n = 7, name = "Set3")[1],"#FFFF00",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)]),c("#5F9EA0","#4169E1","#CD5C5C","#B22222"),
              c("#98F5FF","#8EE5EE","#7AC5CD","#436EEE","#3A5FCD","#27408B","#FF6A6A","#EE6363","#CD5555","#EE2C2C","#CD2626","#8B1A1A"),
              gg_color_hue(2),c("#EC672A","#E5E5E5"),gg_color_hue(20),gg_color_hue(20),
              rep(c(brewer.pal(n = 7, name = "Set3")[1],"#FFFF00",brewer.pal(n = 7, name = "Set3")[c(3,4,7,5,6)]),each=2))
Metadata_column<-c(rep("Annotation",7),rep("Sample",4),rep("Replicate",12),rep("Spleen_localisation",2),rep("Genotype",2),
                   rep("Numbered_SCT_clusters",20),rep("Annotation_detailed",20),rep("Annotation_split_by_genotype",14))
Info_Kevin<-as.data.frame(cbind(All,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results_merge_non_harmony/Annotation/Info_Kevin_",sampleName,".xlsx"))

########################

##### Save object
saveRDS(seuratObj_filtered_diet, file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

##### Read object
seuratObj_filtered_diet<-readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

