## Script for analyzing Xbp1/Ire1KO cDC1 CITE-seq project data
## Muscat DS analysis between WT and DKO cDC1 data

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('ggplot2')
library('openxlsx')
library('cowplot')
library('muscat')
library('purrr')
library('tidyverse')
library('Matrix')
library("limma")
library("RColorBrewer")

library('clusterProfiler')
library('org.Mm.eg.db')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Final/PROCESSED_DATA/")

sampleName<-"SAM_merge"
sampleFolder<-paste0(sampleName,"/")


########################################
##### Some variables
########################################

listLabels<-list(c('SAM05',"SAM06",'SAM07',"SAM08"))

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

### Modified function for robust edgeR
source('~/VIB/DATA/Roos/Cara/pbDSrob.R')
source('~/VIB/DATA/Roos/Cara/pbDS_DESeq2.R')

##### Read final object
seuratObj <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/seuratObj_",sampleName,"_2023.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/diagnostics_",sampleName,"_2023.rds"))

########################################################################################################################

DimPlot(seuratObj, reduction = "SCT_umap", group.by = "annotation_paper_2023")
DimPlot(seuratObj, reduction = "SCT_umap", group.by = "annotation_paper_2023", split.by = "Genotype")

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj@meta.data[["replicate_info"]]

## Combine (same mice for RP and WP)
seuratObj@meta.data[["replicate_muscat"]]<-seuratObj@meta.data[["replicate_info"]]
seuratObj$replicate_muscat<-gsub("SAM05", "WT",seuratObj$replicate_muscat)
seuratObj$replicate_muscat<-gsub("SAM06", "WT",seuratObj$replicate_muscat)
seuratObj$replicate_muscat<-gsub("SAM07", "DKO",seuratObj$replicate_muscat)
seuratObj$replicate_muscat<-gsub("SAM08", "DKO",seuratObj$replicate_muscat)

seuratObj@meta.data$replicate_muscat<-as.factor(seuratObj@meta.data$replicate_muscat)
seuratObj@meta.data$replicate_muscat<-factor(seuratObj@meta.data$replicate_muscat,
                                             levels(seuratObj@meta.data$replicate_muscat)[c(4,5,6,1,2,3)]) #reorder levels

########################################
##### Load data (already filtered) 
########################################

#Convert sparse matrix to matrix 
seuratObj[['SCT']]@counts <- as.matrix(seuratObj[['SCT']]@counts) 
seuratObj[['SCT']]@data <- as.matrix(seuratObj[['SCT']]@data) 

### Convert to sce ### RNA or SCT!!!!!!
sce <- as.SingleCellExperiment(seuratObj, assay = "SCT")

########################################
##### Prepare sce object
########################################

### Create test metadata table
sample_id<-levels(as.factor(seuratObj@meta.data$replicate_muscat))
group_id<-c(rep("WT",3),rep("DKO",3)) 
n_cells<-rep(1,6)

### Create metaData matrix
metaData<-data.frame(sample_id,group_id,n_cells)

### Look through samples (1->8) for amount of cells after filtering
for(i in 1:length(sample_id)){
  toSearch<-sample_id[i]
  metaData[i, "n_cells"]<-length(grep(paste0("\\<",toSearch,"\\>"),seuratObj@meta.data$replicate_muscat)) ##Need to do specific grep. Otherwise fault between Hashtag 1 and Hashtag 10!!!
}

### Retrieve the sizes of the two groups
WT_size<-sum(metaData$n_cells[1:3]) #10649
DKO_size<-sum(metaData$n_cells[4:6]) #12377

### Add group_id to sce object
sce$group_id<-c(rep("WT",WT_size),rep("DKO",DKO_size))

### Add new sample_id to sce object (combination of orig.ident and group_id)
sce$sample_id<-sce$replicate_muscat

### Finish prep object
(sce <- prepSCE(sce, 
                cluster_id = "annotation_paper_2023", # cell population assignments
                group_id = "group_id",   # group IDs (resp/Non_Resp)
                sample_id = "sample_id",    # sample IDs (group+orig.ident)
                drop = TRUE))        # drop all other colData columns

### Store ids separately for easy access
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

############################
##### Visualization
############################

### Calculate number of cells per cluster-sample (can already exclude here or later in aggregate step)
t(table(sce$cluster_id, sce$sample_id))


##############################################################################
##### Differential state analysis (only on sample level for this experiment)
##############################################################################
# Aggregation-based method that acts on *pseudobulk* data. Each gene is tested for state changes in each cluster.
# Thus, a total of #(genes)\times\#(clusters) tests will be performed per comparison of interest.

# First aggregate the measurements for each sample (in each cluster) to obtain pseudobulk data.
# In general, `aggregateData()` will aggregate the data by the `colData` variables specified with argument `by`, 
# and return a `SingleCellExperiment` containing pseudobulk data.  

# For DS analysis, measurements must be aggregated at the cluster-sample level (default `by = c("cluster_id", "sample_id"`). 
# In this case, the returned `SingleCellExperiment` will contain one assay per cluster, where rows = genes and columns = samples. 
# Arguments `assay` and `fun` specify the input data and summary statistic, respectively, to use for aggregation.  
# Default it's applied to the sum of the raw counts, but there are also other possiblities (input data (raw/(log-)normalized counts, 
# CPM ect.) and summary statistics (sum, mean, median))

### Aggregation
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
# pseudobulks for 1st subpopulation
t(head(assay(pb)))

### Pseudobulk-level MDS plot
# Multi-dimensional scaling (MDS) plot of aggregated signal to explore overall sample similarities.
# Ideally, such a representation of the data should separate both clusters and groups from one another. -> Not the case!!
(pb_mds <- pbMDS(pb) + scale_shape_manual(values=1:nlevels(pb$group_id))) #Add function for more than 6 shapes (9 groups)

# Removing cluster-sample instance(s) ‘Ly6c2+ Mature NK cells’-‘DKO_Hashtag3’, ‘Ifitm+ Mature NK cells’-‘WT_Hashtag3’, 
# ‘Ifitm+ Mature NK cells’-‘DKO_Hashtag3’, ‘NK cells’-‘DKO_Hashtag3’, ‘Spp1+Ifng+Ccl3+Ccl4+ Immature NK cells’-‘WT_Hashtag3’, 
# ‘Spp1+Ifng+Ccl3+Ccl4+ Immature NK cells’-‘DKO_Hashtag3’, ‘Gzmc+Ifitm+ Mature NK cells’-‘WT_Hashtag3’, 
# ‘Gzmc+Ifitm+ Mature NK cells’-‘DKO_Hashtag3’, ‘Klra5+ Mature NK cells’-‘DKO_Hashtag3’

dir.create(paste0(sampleFolder,"results_merge_non_harmony/Muscat/"))

Comparison<-"DS_SCT_paper_2023"

png(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/meanVariancePlot_subset_",Comparison,".png"),width = 3500, height = 2000, res = 300)
pb_mds + scale_color_manual(values = RColorBrewer::brewer.pal(7, "Set3"))
dev.off()

### Experimental design
# We can provide `pbDS` with a design matrix capturing the experimental design using `model.matrix` (package `r Rpackage("stats")`), 
# and a contrast matrix that specifies our comparison of interesting using `makeContrasts` from the `r Biocpkg("limma")` package. 
# Alternatively, the comparison(s) of interest (or a list thereof) can be specified with via `coefs` (see `?glmQLFTest` for details). 

# Here, we want to carry out a single comparison of DKO against WT samples, 
# thus placing `"WT"` on the right-hand side as the reference condition.
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id) #model.matrix(~ 0 + ei$group_id + ei$gender) Include covariates here!!
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("DKO-WT",
                          levels = mm)

# run DS analysis
### DeSeq2 method ###
res_DESeq2<-pbDS(pb, design = mm, contrast = contrast,method = "DESeq2")

## Save res object
saveRDS(res_DESeq2,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds"))

## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

## Create loop for all comparisons
for (i in 1:length(res$table)){
  tbl <- res$table[[i]]
  Comp<-names(res$table)[i]
  write.xlsx(tbl, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_full_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })
  
  # filter clint
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.5)
    dplyr::arrange(u, p_adj.loc)
  })
  
  saveRDS(tbl_fil,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
  saveRDS(tbl_fil_clint,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_clint, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

## Relax the logFC cutoff further (final version for the paper??)
for (i in 1:length(res$table)){
  tbl <- res$table[[i]]
  Comp<-names(res$table)[i]

  # filter clint
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.4)
    dplyr::arrange(u, p_adj.loc)
  })

  saveRDS(tbl_fil_clint,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/tbl_fil_muscat_DESeq2_relaxed_v2_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_clint, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_fil_muscat_DESeq2_relaxed_v2_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

## Relax the logFC cutoff further + split up by logFC (UP vs DOWN) (final version 04/2024 for the paper!!!)
for (i in 1:length(res$table)){
  tbl <- res$table[[i]]
  Comp<-names(res$table)[i]
  
  # filter clint
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, logFC > 0.4)
    dplyr::arrange(u, p_adj.loc)
  })
  
  tbl_fil_clint2 <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, logFC < -0.4)
    dplyr::arrange(u, p_adj.loc)
  })
  
  names(tbl_fil_clint)<-paste0(names(tbl_fil_clint),"_UP")
  names(tbl_fil_clint2)<-paste0(names(tbl_fil_clint2),"_DOWN")
  tbl_fil_final<-c(tbl_fil_clint,tbl_fil_clint2)
  tbl_fil_final<-tbl_fil_final[c(8,1,9,2,10,3,11,4,12,5,13,6,14,7)]
  
  saveRDS(tbl_fil_final,file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/tbl_fil_muscat_DESeq2_relaxed_v2_split_paper_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_final, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_fil_muscat_DESeq2_relaxed_v2_split_paper_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

### Extra ###
tbl_freq<- resDS(sce, res, frq = TRUE)
write.xlsx(tbl_freq, paste0(sampleFolder,"results_merge_non_harmony/Muscat/tbl_full_freq_muscat_DESeq2_new_",Comparison,"_",sampleName,".xlsx"))

## Get list of cell populations
listCells<-levels(sce$cluster_id)[c(1:7)] #Empty lists!!

names(tbl_fil_clint)

pbHeatmap(sce, res, top_n = 20)

for (i in 1:length(res$table)){
  H2 <- pbHeatmap(sce, res, top_n = 25, c = names(res$table[i]), k = listCells, normalize = T) 
  
  Comp<-names(res$table)[i]
  
  pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Heatmap_Overview_DESeq2_new_normalized_",Comp,".pdf"), width = 8, height = 15)
  print(H2)
  dev.off()
  
}

################
################
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Muscat/Violin_plots/"))
dir.create(paste0(sampleFolder,"results_merge_non_harmony/Muscat/Heatmaps/"))

##Regular loop:
for (i in 1:length(listCells)) {
  gene_amount<-length(tbl_fil[[listCells[i]]][,"gene"])
  if (gene_amount > 6) {
    gene_amount<-6
  }
  p <- plotExpression(sce[, sce$cluster_id == listCells[i]],
                      features = tbl_fil[[listCells[i]]][seq_len(gene_amount),"gene"],
                      x = "sample_id", colour_by = "group_id", ncol = 3) +
    guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(p, filename=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Violin_plots/Violin_plot_new_",listCells[i],"_DESeq2_",Comparison,".png")) # width = 1200, height = 800
  
  gene_amount<-length(tbl_fil[[listCells[i]]][,"gene"])
  if (gene_amount > 20) {
    gene_amount<-20
  }
  
  H<- pbHeatmap(sce, res, g = tbl_fil[[listCells[i]]][seq_len(gene_amount),"gene"] , k = listCells[i], normalize = T, lfc = 1) #Add lfc
  
  png(filename=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Heatmaps/Heatmap_new_",listCells[i],"_DESeq2_normalized_",Comparison,".png"), width=1500, height=1000, res = 300)
  print(H)
  dev.off()
}

# Paper loop (Added col argument)
for (i in 1:length(listCells)) {
  gene_amount<-length(tbl_fil[[listCells[i]]][,"gene"])
  if (gene_amount > 50) {
    gene_amount<-50
  }
  
  H<- pbHeatmap(sce, res, g = tbl_fil[[listCells[i]]][seq_len(gene_amount),"gene"] , top_n = 50, k = listCells[i], sort_by = "logFC", normalize = T, lfc = 1,
                col = c(rev(c('#d1e5f0','#67a9cf','#2166ac')),"white",rev(c('#b2182b','#ef8a62','#fddbc7'))))
  
  pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Heatmaps/Heatmap_logFC_",listCells[i],"_DESeq2_normalized_top50_",Comparison,".pdf"), width=5, height=6)
  print(H)
  dev.off()
}

# Paper loop final version (Added col argument)
for (i in 1:length(listCells)) {
  gene_amount<-length(tbl_fil[[listCells[i]]][,"gene"])
  if (gene_amount > 50) {
    gene_amount<-50
  }
  
  H<- pbHeatmap(sce, res, g = tbl_fil[[listCells[i]]][seq_len(gene_amount),"gene"] , top_n = 50, k = listCells[i], sort_by = "p_adj.loc", normalize = T, lfc = 1,
                col = c(rev(c('#d1e5f0','#67a9cf','#2166ac')),"white",rev(c('#b2182b','#ef8a62','#fddbc7'))))
  
  pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Heatmaps/Heatmap_adjPV_",listCells[i],"_DESeq2_normalized_top50_",Comparison,"_",sampleName,".pdf"), width=5, height=6)
  print(H)
  dev.off()
}


######################################

## Extra 2025: create specific heatmap for paper!!

## EM cDC1s + LM cDC1s
## Subset average seuratobject 
seuratObj_average_subset<-subset(seuratObj_average, idents = levels(Idents(seuratObj_average))[c(31:42)])

## Subset DEG list without Rps/Rpl genes
Paper_genelist_Mature<-c("Ggta1","Nabp1","Bcl2l11","Nfkbia","Bcl2l14","Ccr7","Cblb","Papss2","Tnfrsf4",
                         "Mxd1","Traf1","Clec2i","Irf1","Marcksl1","Asprv1","Snn","Rras2","Rasip1","Tmem176a",
                         "Dennd3","Tmem39a","Bcl3","Il2rg","Il4ra","Slc2a6","Ccl22","Ldlr","Kcnk6","Lamp1","Gbp8",
                         "Zfp36l1","St8sia6","Chka","Slco5a1","Rasal2","Mical3","Mmp23","Eno3","Sqle","Cyp51","Creld1",
                         "Ccr7","Fscn1","Tmem176a","Tmem176b","Slco5a1","Nudt17","Apol7c","Apoe",
                         "Apol10b","Abcg1","Il4i1","Mical3","H2-M2","Dnase1l3")

## Heatmap mature cDC1s
# Colset<-c(rep('royalblue',3),
#           rep('firebrick',3))
Colset<-rep(brewer.pal(n = 6, name = "Set3"), each = 6)
Colset_Mature<-Colset[-c(1:24)] 
H1 <- DoHeatmap(seuratObj_average_subset, assay = "SCT", #cells = levels(Idents(seuratObj_average))[c(1:6)],
                features = Paper_genelist_Mature, size = 3,
                draw.lines = FALSE, group.colors = Colset_Mature,disp.min = -2.5, disp.max = 2.5) + 
  # scale_color_manual(values = Colset_EM) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Seurat_heatmaps/Average_expression_heatmap_Seurat_paper_DESeq2_Maturation_genes_EM_and_LM_",sampleName,".pdf"), width = 12, height = 15)
print(H1)
dev.off()

## Heatmap all cDC1s
Colset<-rep(brewer.pal(n = 7, name = "Set3"), each = 6)
Colset_All<-Colset[c(1:24,37:42,25:36)] 
H1 <- DoHeatmap(seuratObj_average, assay = "SCT", #cells = levels(Idents(seuratObj_average))[c(1:6)],
                features = Paper_genelist_Mature, size = 3,
                draw.lines = FALSE, group.colors = Colset_All,disp.min = -2.5, disp.max = 2.5) + 
  # scale_color_manual(values = Colset_EM) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Seurat_heatmaps/Average_expression_heatmap_Seurat_paper_DESeq2_Maturation_genes_all_cDC1s_",sampleName,".pdf"), width = 20, height = 15)
print(H1)
dev.off()

## ComplexHeatmap
library(ComplexHeatmap)
features <- unique(Paper_genelist_Mature)

## Update genelist: split in DEG and non-DEG!!!
DEG_features<-unique(c(intersect(features,tbl_fil_clint_relaxed_v2$`Late mature cDC1s`$gene),
                       intersect(features,tbl_fil_clint_relaxed_v2$`Early mature cDC1s`$gene))) 
nonDEG_features<-setdiff(features,DEG_features)

# Extract the expression data for the selected genes from the scaled layer
expr_matrix <- as.matrix(
  GetAssayData(
    seuratObj_average, 
    slot = "scale.data", 
    assay = "SCT"
  )[
    c(nonDEG_features,DEG_features), #features
  ]
)

# Prepare cell type annotations
cell_types <- Idents(seuratObj_average)
# Create vector with unique and sorted cell types
unique_cell_types <- levels(cell_types)

## Add genotype info
seuratObj_average@meta.data$Genotype<-c(rep("WT",3),rep("DKO",3))
seuratObj_average@meta.data$Genotype<-factor(seuratObj_average@meta.data$Genotype, levels = c("WT","DKO"))

# Create color palette for cell types
cell_type_colors <- setNames(
  Colset_All,
  unique_cell_types
)

# Create color palette for cell types
Genotype_colors <- setNames(
  rep(c(rep("#E5E5E5",3),rep("#EC672A",3)),7),
  seuratObj_average@meta.data$Genotype
)

# Recreate the Seurat expression scale
palette_expression_level = circlize::colorRamp2(
  c(min(expr_matrix), median(expr_matrix), min(2.5, max(expr_matrix))),
  c("#2166ac", "white", '#b2182b'))

# Create the heatmap with ComplexHeatmap
ht = Heatmap(
  matrix = expr_matrix,
  row_order = c(nonDEG_features,DEG_features), #features,
  name = "expression",
  column_split = factor(cell_types, levels = unique_cell_types),
  row_split = c(rep("1.non-DEG",30), rep("2. DEG",21)),
  # Do not cluster rows or columns otherwise row_order/column_order are ignored
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_row_dend = FALSE,
  cluster_column_slices = TRUE,
  column_gap = unit(0.5, "mm"),
  # column_names_rot = 45,
  # column_names_side = "top",
  column_names_max_height = unit(5, "cm"),
  row_names_gp = gpar(fontsize = 5),
  # The column title can be either a title for all the columns
  # or a vector with a name for every column
  column_title = unique_cell_types,
  column_title_gp = gpar(fontsize = 2),
  column_title_rot = 90,
  # Add the annotation bars on top
  top_annotation = HeatmapAnnotation(
    # Give a name to the annotation, not used on the actual plot
    name = "Clusters",
    # Annotate the columns
    which = "column",
    # First annotation - clusters
    cluster = cell_types,
    # Second annotation - patients
    genotype = seuratObj_average@meta.data$Genotype,
    # Color palettes
    col = list(cluster = cell_type_colors, genotype = Genotype_colors),
    show_legend = c(FALSE, TRUE),
    gap = unit(1, "mm"),
    simple_anno_size = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 6)
  ),
  # Color scale for expression
  col = palette_expression_level,
  # Otherwise you get random vertical white lines where there should not be any
  use_raster = FALSE
)
# Put the legends on the right in a single column
pdf(file=paste0(sampleFolder,"results_merge_non_harmony/Muscat/Seurat_heatmaps/Complexheatmap_paper_DESeq2_Maturation_genes_all_cDC1s_updated",sampleName,".pdf"), width = 10, height = 8)
draw(ht, merge_legend = TRUE)
dev.off()

## Compare results bulk and CITEseq logFC
Bulk_results<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/triwiseResults/summary_allgenes_clint.xlsx", sheet = "Ire1KO_vs_WT")
Pseudobulk_DEG_results<-tbl_fil_clint_relaxed_v2$`Late mature cDC1s`
Pseudobulk_results<-res$table$`DKO-WT`$`Late mature cDC1s`

rownames(Bulk_results)<-Bulk_results$gene
rownames(Pseudobulk_results)<-Pseudobulk_results$gene

Bulk_results[c(nonDEG_features,DEG_features),c(7,1,2,5)]
Pseudobulk_results[c(nonDEG_features,DEG_features),c(2,3,4,8)]

write.xlsx(cbind(Bulk_results[c(nonDEG_features,DEG_features),c(7,1,2,5)],
                 Pseudobulk_results[c(nonDEG_features,DEG_features),c(2,3,4,8)]),
           paste0(sampleFolder,"results_merge_non_harmony/Muscat/Seurat_heatmaps/Complexheatmap_paper_DESeq2_Maturation_genes_all_cDC1s_updated",sampleName,".xlsx"))

########################################################################################################
########################################################################################################

#################
### GSEA in R ###
#################

#Create dir
dir.create(paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/"))

#Load markers
## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

for (i in 1:length(res$table)){
  
  tbl <- res$table[[i]]
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })
  
  Comp<-names(res$table)[i]
  
  ## All cell clusters
  EnrichGO<-tibble::lst()
  for (cell_pop in 1:length(tbl_fil)) { 
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #EnrichGO:
    EnrichGO[cell_pop]<-enrichGO(
      as.character(tbl_fil[[cell_pop]]$gene),
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO, paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_results_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  for (cell_pop in (1:length(tbl_fil))) { 
    ## Create dotplot top 10 each category
    D_background<-dotplot(EnrichGO[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    
    ## Save files
    pdf(file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/Dotplot_background_scRNAseq_",names(tbl_fil)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
    print(D_background)
    dev.off()
    
    write.xlsx(EnrichGO[[cell_pop]]@result,file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_background_scRNAseq_",names(tbl_fil)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
  
}

###################

## Relaxed version 1

#Load markers
## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

for (i in 1:length(res$table)){
  
  tbl <- res$table[[i]]
  
  # filter clint
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.5)
    dplyr::arrange(u, p_adj.loc)
  })
  
  Comp<-names(res$table)[i]
  
  ## All cell clusters
  EnrichGO<-tibble::lst()
  for (cell_pop in 1:length(tbl_fil_clint)) { 
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #EnrichGO:
    EnrichGO[cell_pop]<-enrichGO(
      as.character(tbl_fil_clint[[cell_pop]]$gene),
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO, paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_relaxed_results_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  ## Read results
  # EnrichGO <- readRDS(paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_results_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  
  for (cell_pop in (1:length(tbl_fil_clint))) { 
    ## Create dotplot top 10 each category
    D_background<-dotplot(EnrichGO[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    
    ## Save files
    pdf(file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/Dotplot_relaxed_results_background_scRNAseq_",names(tbl_fil_clint)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
    print(D_background)
    dev.off()
    
    write.xlsx(EnrichGO[[cell_pop]]@result,file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_relaxed_results_background_scRNAseq_",names(tbl_fil_clint)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
  
}

###################

## Relaxed version 2

#Load markers
## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results_merge_non_harmony/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

for (i in 1:length(res$table)){
  
  tbl <- res$table[[i]]
  
  # filter clint
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.4)
    dplyr::arrange(u, p_adj.loc)
  })
  
  Comp<-names(res$table)[i]
  
  ## All cell clusters
  EnrichGO<-tibble::lst()
  for (cell_pop in 1:length(tbl_fil_clint)) { 
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #Genelist
    Genelist_scRNAseq<-as.character(tbl_fil_clint[[cell_pop]]$gene)
    #EnrichGO:
    EnrichGO[cell_pop]<-enrichGO(
      Genelist_scRNAseq,
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO, paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_relaxed_results_v2_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  for (cell_pop in (1:length(tbl_fil_clint))) { 
    ## Create dotplot top 10 each category
    D_background<-dotplot(EnrichGO[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    
    ## Save files
    pdf(file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/Dotplot_relaxed_results_v2_background_scRNAseq_",names(tbl_fil_clint)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
    print(D_background)
    dev.off()
    
    write.xlsx(EnrichGO[[cell_pop]]@result,file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_relaxed_results_v2_background_scRNAseq_",names(tbl_fil_clint)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
}

## Repeat relaxed version 2: only look at early and late mature cDC1s and split up by logFC (January 2024)
for (i in 1:length(res$table)){
  
  tbl <- res$table[[i]]
  
  # filter clint
  tbl_fil_clint_UP <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, logFC > 0.4)
    dplyr::arrange(u, p_adj.loc)
  })
  
  tbl_fil_clint_DOWN <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, logFC < -0.4)
    dplyr::arrange(u, p_adj.loc)
  })
  
  Comp<-names(res$table)[i]
  
  ## UP cell clusters
  EnrichGO_UP<-tibble::lst()
  for (cell_pop in 5:length(tbl_fil_clint_UP)) { 
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #Genelist
    Genelist_scRNAseq<-as.character(tbl_fil_clint_UP[[cell_pop]]$gene)
    #EnrichGO:
    EnrichGO_UP[cell_pop]<-enrichGO(
      Genelist_scRNAseq,
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO_UP, paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_UP_relaxed_results_v2_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  for (cell_pop in (5:length(tbl_fil_clint_UP))) { 
    ## Create dotplot top 10 each category
    D_background<-dotplot(EnrichGO_UP[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    
    ## Save files
    pdf(file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/Dotplot_relaxed_results_v2_UP_background_scRNAseq_",names(tbl_fil_clint_UP)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
    print(D_background)
    dev.off()
    
    write.xlsx(EnrichGO_UP[[cell_pop]]@result,file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_UP_relaxed_results_v2_background_scRNAseq_",names(tbl_fil_clint_UP)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
  
  ## DOWN cell clusters
  EnrichGO_DOWN<-tibble::lst()
  for (cell_pop in 5:length(tbl_fil_clint_DOWN)) { 
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #Genelist
    Genelist_scRNAseq<-as.character(tbl_fil_clint_DOWN[[cell_pop]]$gene)
    #EnrichGO:
    EnrichGO_DOWN[cell_pop]<-enrichGO(
      Genelist_scRNAseq,
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO_DOWN, paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_DOWN_relaxed_results_v2_background_scRNAseq_",Comp,"_",sampleName,".rds"))
  
  for (cell_pop in (5:length(tbl_fil_clint_DOWN))) { 
    if (nrow(EnrichGO_DOWN[[cell_pop]]) > 1){ #Empty list!!
      ## Create dotplot top 10 each category
      D_background<-dotplot(EnrichGO_DOWN[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
      
      ## Save files
      pdf(file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/Dotplot_relaxed_results_v2_DOWN_background_scRNAseq_",names(tbl_fil_clint_DOWN)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
      print(D_background)
      dev.off()
    }
    
    write.xlsx(EnrichGO_DOWN[[cell_pop]]@result,file=paste0(sampleFolder,"results_merge_non_harmony/GSEA_Muscat/EnrichGO_DOWN_relaxed_results_v2_background_scRNAseq_",names(tbl_fil_clint_DOWN)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
  
}