## Script for analyzing Xbp1/Ire1KO cDC1 bulk RNA-seq project data
## limma-edgeR DEG pipeline combined with triwise and GO analysis

library("limma")
library("edgeR")
library("ggplot2")
library("triwise")
library("htmlwidgets")
library("rgl")
library("RColorBrewer")
library("gplots")
library('xlsx')
library('openxlsx')

library("org.Mm.eg.db")
library("GO.db")
library("dplyr")
library('gridExtra')

library('readxl')
library('stringr')

library('pheatmap')
library('grid')
library('gplots')
library('ggrepel')

########################################
##### Functions
########################################

###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###First letter upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}

################################################################################
######### LOAD META DATA
################################################################################

getwd()
setwd("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/scripts")
savePath <- "/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/"

### Load raw counts
countData<-read.table(file="../rawData/outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE) 

### Load meta data
colData<-read.table(file="../rawData/outputHTSeqCount/outputFilesHTSeqCount_Ire1_Xbp1_KOmice/metadata_v2.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE) 
rownames(colData)<-colData$fileName
colData<-colData[,-1]

### Only cDC1, without 'Xbp1_Ire1_WT_cDC1_1' and 'Xbp1_KO_cDC1_2' and 'Xbp1_Ire1_KO_cDC1_4'
colData<-colData[c(2:7,9,11:12),]
countData<-countData[,rownames(colData)] 

################################################################################
########## CREATE OBJECT
################################################################################

y <- DGEList(counts = countData)

################################################################################
########## FILTER DATA
################################################################################

##### Filter low count genes
## always work with count-per-million (CPM) instead of raw counts
## Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library
## Imagine that the lowest lib size is around 6 million reads => threshold is set on CPM>1
## But what if lib size is 20 million? Here CPM of 1 means 20 counts. Do we still use 5 counts and then the CPM cut off
## will be 0.25 or is the threshold of 5 counts a cut off for lib sizes of around 5-6 million? Then we need to put the
## cut off on 20 for lib sizes around 20 million and so use a CPM of 1.
## Threshold needs to be true in at least x samples. x is always the lowest number of replicates.
## for example: 3 samples with each 2 replicates => x set on 2
## => This ensures that a gene will be retained if it is only expressed in both replicates of a certain group

## Do filtering
yNoFilter<-y
myCpm<-cpm(y)

keep = rowSums(cpm(y)>1) >= 3
y = y[keep,]

##### Reset lib sizes
y$samples$lib.size = colSums(y$counts)
y$samples

################################################################################
########## NORMALISATION
################################################################################

##### Scale normalisation
yNoNormalisation<-y
y <- calcNormFactors(y)

##### MDS-plot
theColors<-c(rep("royalblue",3),rep("darkgreen",3),rep("limegreen",3))

plotMDS(y, cex=0.8, col=theColors)

################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
TS <- paste(colData$cell, colData$condition, sep=".")
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
v <- voom(y, design, plot = TRUE)

expTable<-v$E
saveRDS(expTable,"../results/expTable.rds")

###########################################
###########################################

## Update 2023 for upload to GEO
dir.create("../results/GEO")

# Check names
all(colnames(countData) == colnames(expTable))

write.table(countData, "../results/GEO/Bulk_RNA_seq_raw_counts_cDC1s_Ire1_Xbp1_experiment.txt", sep="\t")
write.table(expTable, "../results/GEO/Bulk_RNA_seq_normalized_counts_cDC1s_Ire1_Xbp1_experiment.txt", sep="\t")

################################################################################
########## CHECK FILTERING AND NORMALISATION
################################################################################

#################### BARPLOT ####################

##### Barplot lib sizes raw counts
par(mar = c(9,3,3,1)) #more margin: bottom, left, top, right
bp<-barplot(yNoFilter$samples$lib.size*1e-6,axisnames=FALSE,main="Barplot lib sizes of raw counts",ylab="Library size (millions)")
axis(1, labels=rownames(yNoFilter$samples), at = bp, las=2, cex.axis=0.7)

##### Barplot lib sizes normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
bp<-barplot(colSums(countData_norm)*1e-6,axisnames=FALSE,main="Barplot lib sizes of normalised counts",ylab="Library size (millions)")
axis(1, labels=colnames(countData_norm), at = bp, las=2, cex.axis=0.7)


#################### BOXPLOT ####################
y2<-y
y2$samples$norm.factors<-1
y2$counts[,1]<-ceiling(y2$counts[,1]*0.05)
y2$counts[,2]<-y2$counts[,2]*5

par(mfrow=c(1,2), mar = c(12,4,3,1)) #more margin: bottom, left, top, right
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Unnormalised data", ylab="log-cpm")

y2<-calcNormFactors(y2)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Normalised data", ylab="log-cpm")


#################### DENSITY PLOT ####################
col <- topo.colors(nrow(colData))

par(mfrow=c(1,2))
### Plot log2-CPM values of each sample before filtering
theCpmNoFilter<-cpm(yNoFilter, log=TRUE)

plot(density(theCpmNoFilter[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="raw data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpmNoFilter)[i]))
  
  den <- density(theCpmNoFilter[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}

### Plot log2-CPM values of each sample after filtering (and normalisation)
theCpm<-cpm(y, log=TRUE)

plot(density(theCpm[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="filtered data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpm)[i]))
  
  den <- density(theCpm[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
par(mfrow=c(1,1))

#################### HISTOGRAM OF EXPTABLE ####################

par(mfrow=c(1,2))
### Histogram
hist(expTable)

### Density plot
plot(density(expTable[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="Density of expTable", xlab="log2")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0("Add line for sample ",colnames(expTable)[i]))
  
  den <- density(expTable[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
par(mfrow=c(1,1))

################################################################################
########## PCA
################################################################################

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab="pc1", ylab="pc2", zlab="pc3")
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.3)

################################################################################
########## CORRELATION HEATMAP SAMPLES
################################################################################

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

##heatmap 1: based on distance
distsRL <- dist(t(expTable),method="euclidean")
hc <- hclust(distsRL,method="ward.D")

heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13), cexRow=0.9, cexCol=0.9)


##heatmap 2: based on correlation
cm=cor(expTable)

distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")


heatmap.2(cm, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=hmcol, margin=c(13, 13), cexRow=0.9, cexCol=0.9)


################################################################################
########## GET DE GENES
################################################################################

#### Fit linear model on data
fit <- lmFit(v, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=cDC1.XBPKO-cDC1.WT, group2=cDC1.IREKO-cDC1.WT, group3=cDC1.IREKO-cDC1.XBPKO, levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)


#### Adjust P-values via Benjamini-Hochberg
##############################
##### Xbp1 KO vs WT
##############################
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.05,1)
dim(DEgenesGroup1)
##418

##############################
##### Ire1 KO vs WT
##############################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,1)
dim(DEgenesGroup2)
##338

##############################
##### Ire1 KO vs Xbp1 KO
##############################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.05,1)
dim(DEgenesGroup3)
##413


### Write results
write.xlsx(DEgenesGroup1, "../results/triwiseResults/summary_DEgenes_new.xls", "Xbp1KO_vs_WT")
write.xlsx(DEgenesGroup2, "../results/triwiseResults/summary_DEgenes_new.xls", "Ire1KO_vs_WT", append=T)
write.xlsx(DEgenesGroup3, "../results/triwiseResults/summary_DEgenes_new.xls", "Ire1KO_vs_Xbp1KO", append=T)

#####Full list (10/06/20)
coef<-c(1:ncol(cont.matrix))
listallgenes<-list()

for(i in 1:length(coef)){
  allGenes<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=coef[i])
  listallgenes[[i]]<-allGenes
}


names(listallgenes)<-c("Xbp1KO_vs_WT","Ire1KO_vs_WT",'Ire1KO_vs_Xbp1KO')
listallgenes<-lapply(listallgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})

write.xlsx(listallgenes, file = "../results/triwiseResults/summary_allgenes_clint.xlsx")

####################################################################################################
############################## ----------- TRIWISE PLOTS ------------ ##############################
####################################################################################################

##### All DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2),rownames(DEgenesGroup3)))
length(allDEgenes)
##790

######################################################################
########## PREPARATION
######################################################################
wantedColors7<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red",
                nodiffSet2="limegreen",diffSet2="darkgreen", nodiffSet3="lightblue", diffSet3="darkblue",
                nodiffSet4="darkorchid1", diffSet4="darkmagenta", nodiffSet5="lightsalmon", diffSet5="darkorange",
                nodiffSet6="lightskyblue", diffSet6="turquoise4", nodiffSet7="bisque", diffSet7="burlywood4")

##### Triwise plots
colsWT<-c(grep("WT",colnames(expTable)), grep("fl_fl",colnames(expTable)))
colsIre1KO<-grep("IRE1_KO",colnames(expTable))
colsXpb1KO<-grep("XBP1_KO",colnames(expTable))

WT_mean<-apply(expTable[,colsWT],1,mean)
Ire1KO_mean<-apply(expTable[,colsIre1KO],1,mean)
Xbp1KO_mean<-apply(expTable[,colsXpb1KO],1,mean)

expTable_mean<-cbind(WT_mean,Ire1KO_mean,Xbp1KO_mean)
colnames(expTable_mean)<-c('WT','Ire1KO','Xbp1KO')

######################################################################
########## TRIWISE PLOTS
######################################################################
theBarycoords<-transformBarycentric(expTable_mean)

barycoords = theBarycoords
Gdiffexp = allDEgenes

saveRDS(barycoords, "../results/barycoords_new.rds")
saveRDS(Gdiffexp, "../results/Gdiffexp_new.rds")

###Triwise plot
p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors7, showlabels = F, sizevalues = 50)
print(p)
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_new.png",dpi=300)

###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_new.png",dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_rosePlot.pdf"))
plotRoseplot(barycoords, Gdiffexp, Goi=lipidGenes, showlabels = F) ##Use genelist before intersect (but not necessary)!
dev.off()

###Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE)
print(p)
saveWidget(p,file="C:/DATA/RNA-seq_Simon/results/triwiseResults/interactiveTriwisePlot_new.html") ##needs full path!

###Seperate the DE genes of the different comparisons
#set1 (red)=DE genes only present in Xbp1 KO vs WT
#set2 (green)=DE genes only present in Ire1 KO vs WT
#set3 (blue)=DE genes only present in Ire1 KO vs Xpb1 KO
#set4 (purple)=DE genes present in all 3 comparisons
#set5 (orange)=DE genes present in Xbp1 KO vs WT and Ire1 KO vs WT
#set6 (turquoise)=DE genes present in Ire1 KO vs WT and Ire1 KO vs Xbp1 KO
#set7 (brown)=DE genes present in Xbp1 KO vs WT and Ire1 KO vs Xbp1 KO

genesSet5tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup2))
genesSet6tmp<-intersect(rownames(DEgenesGroup2),rownames(DEgenesGroup3))
genesSet7tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup3))
genesSet4<-intersect(genesSet5tmp, rownames(DEgenesGroup3))
genesSet5<-setdiff(genesSet5tmp,genesSet4)
genesSet6<-setdiff(genesSet6tmp,genesSet4)
genesSet7<-setdiff(genesSet7tmp,genesSet4)
genesSet1<-setdiff(rownames(DEgenesGroup1),unique(c(rownames(DEgenesGroup2),rownames(DEgenesGroup3))))
genesSet2<-setdiff(rownames(DEgenesGroup2),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup3))))
genesSet3<-setdiff(rownames(DEgenesGroup3),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2))))

genesSet1matrix<-DEgenesGroup1[genesSet1,c(1,2,5)]
genesSet2matrix<-DEgenesGroup2[genesSet2,c(1,2,5)]
genesSet3matrix<-DEgenesGroup3[genesSet3,c(1,2,5)]
genesSet4matrix<-cbind(DEgenesGroup1[genesSet4,c(1,2,5)],DEgenesGroup2[genesSet4,c(1,2,5)],DEgenesGroup3[genesSet4,c(1,2,5)],
                       'meanLogFC'=apply(cbind(DEgenesGroup1[genesSet4,1],DEgenesGroup2[genesSet4,1],DEgenesGroup3[genesSet4,1]),1,mean))
genesSet5matrix<-cbind(DEgenesGroup1[genesSet5,c(1,2,5)],DEgenesGroup2[genesSet5,c(1,2,5)],
                       'meanLogFC'=apply(cbind(DEgenesGroup1[genesSet5,1],DEgenesGroup2[genesSet5,1]),1,mean))
genesSet6matrix<-cbind(DEgenesGroup2[genesSet6,c(1,2,5)],DEgenesGroup3[genesSet6,c(1,2,5)],
                       'meanLogFC'=apply(cbind(DEgenesGroup2[genesSet6,1],DEgenesGroup3[genesSet6,1]),1,mean))
genesSet7matrix<-cbind(DEgenesGroup1[genesSet7,c(1,2,5)],DEgenesGroup3[genesSet7,c(1,2,5)],
                       'meanLogFC'=apply(cbind(DEgenesGroup1[genesSet7,1],DEgenesGroup3[genesSet7,1]*-1),1,mean)) #Need to multiply with -1 to correct for Xbp1 KO switching sides of comparison.
### Write results
write.xlsx(genesSet1matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "red")
write.xlsx(genesSet2matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "green", append=T)
write.xlsx(genesSet3matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "blue", append=T)
write.xlsx(genesSet4matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "violet", append=T)
write.xlsx(genesSet5matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "gold", append=T)
write.xlsx(genesSet6matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "cyan", append=T)
write.xlsx(genesSet7matrix, "../results/triwiseResults/triwisePlot_withColors_newcolors.xls", "purple", append=T)

p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F, Goi=list(Set1=genesSet1,Set2=genesSet2,Set3=genesSet3,Set4=genesSet4,
                                                                            Set5=genesSet5,Set6=genesSet6,Set7=genesSet7), colorvalues=wantedColors7) +
  theme(legend.position = "none")
print(p)
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_withColors_newcolors.png",dpi=300)

myList<-list("set1"=genesSet1,"set2"=genesSet2,"set3"=genesSet3,"set4"=genesSet4,"set5"=genesSet5,"set6"=genesSet6,"set7"=genesSet7)
saveRDS(myList,file="../results/listGenes_coloredTriwise_new.rds")

######################################################################
########## FIND GO TERMS
######################################################################

# http://saeyslab.github.io/triwise/vignette.html

### Create list with GO terms
gsets = AnnotationDbi::as.list(org.Mm.egGO2ALLEGS)
symbConv<-unlist(as.list(org.Mm.egSYMBOL))
symbConv2<-names(symbConv)
names(symbConv2)<-symbConv
allGenes<-symbConv2[rownames(expTable_mean)]
allGenes<-allGenes[-1]

##remove all the genes not present in expTable
gsets = sapply(gsets, function(gset) intersect(allGenes, unique(as.character(gset))))
##only keep the BP GO terms
gsetindex = dplyr::bind_rows(lapply(AnnotationDbi::as.list(GOTERM[names(gsets)]), function(goinfo) {
  tibble(name=Term(goinfo), definition=Definition(goinfo), ontology=Ontology(goinfo), gsetid = GOID(goinfo))
  }))
gsets = gsets[gsetindex %>% filter(ontology == "BP") %>% .$gsetid]
##convert EntrezIDs into gene symbols
gsetsRes<-lapply(gsets, function(gset){ gset<-symbConv[gset] })
gsets<-gsetsRes
saveRDS(gsets, "../results/gsets_new.rds")
saveRDS(gsetindex, "../results/gsetindex_new.rds")

# gsets<-readRDS("../results/gsets_new.rds") 
# gsetindex<-readRDS("../results/gsetindex_new.rds")


### Test if GO terms are upregulated in a certain direction (using rank statistic)
scoresFull = testUnidirectionality(barycoords, gsets, Gdiffexp, statistic = "rank", nsamples=1e+6) #mc.cores = 8 : removed this, not supported in Windows.
saveRDS(scoresFull,file="../results/scores_GOterms_newgset_new.rds")
# scoresFull <- readRDS("../results/scores_GOterms_newgset_new.rds")

### Test if GO terms are upregulated in a certain direction (using r statistic)
scoresFull = testUnidirectionality(barycoords, gsets, Gdiffexp, statistic = "r", nsamples=1e+6) #mc.cores = 8 : removed this, not supported in Windows.
saveRDS(scoresFull,file="../results/scores_GOterms_newgset_r_statistic.rds")
# scoresFull <- readRDS("../results/scores_GOterms_newgset_r_statistic.rds")

### Test if GO terms are upregulated in a certain direction (using diffexp statistic)
scoresFull = testUnidirectionality(barycoords, gsets, Gdiffexp, statistic = "diffexp", nsamples=1e+6) #mc.cores = 8 : removed this, not supported in Windows.
saveRDS(scoresFull,file="../results/scores_GOterms_newgset_diffexp_statistic.rds")
# scoresFull <- readRDS("../results/scores_GOterms_newgset_diffexp_statistic.rds")

scores = left_join(scoresFull, gsetindex, by="gsetid") %>% filter(qval < 0.05) %>% arrange(qval, z)
dim(scores)

### Remove the redundancy of GO results
scores = scores[(scores$qval < 0.05) & (scores$z > 3500), ] #With rank statistic
scores = scores[(scores$qval < 0.05) & (scores$z > 0.15), ] #With r statistic
scores = scores[(scores$qval < 0.05) & (scores$z > 0.05), ] #With diffexp statistic

scores$redundancy = estimateRedundancy(scores, gsets, Gdiffexp)
dim(scores)
## 102 with rank and 3500 cutoff
## 115 with r and 0.15 cutoff
## 112 with diffexp and 0.05 cutoff

### Plot the top 20 found GO terms
plotPvalplot(scores %>% top_n(20, redundancy), colnames(expTable_mean))

### Plot all GO terms
plotPvalplot(scores, colnames(expTable_mean))

scoresPart1<-filter(scores, angle >= 0, angle < 1) %>% .[order(.$qval),]
scoresPart2<-filter(scores, angle >= 1, angle < 2) %>% .[order(.$qval),]
scoresPart3<-filter(scores, angle >= 2, angle < 3) %>% .[order(.$qval),]
scoresPart4<-filter(scores, angle >= 3, angle < 4) %>% .[order(.$qval),]
scoresPart5<-filter(scores, angle >= 4, angle < 5) %>% .[order(.$qval),]
scoresPart6<-filter(scores, angle >= 5, angle < 6) %>% .[order(.$qval),]
scoresPart7<-filter(scores, angle >= 6, angle < 7) %>% .[order(.$qval),]
myList<-list("angle_0-1"=scoresPart1,"angle_1-2"=scoresPart2,"angle_2-3"=scoresPart3,"angle_3-4"=scoresPart4,"angle_4-5"=scoresPart5,
             "angle_5-6"=scoresPart6,"angle_6-7"=scoresPart7)
scorePlots<-tibble::lst(scoresPart1,scoresPart2,scoresPart3,scoresPart4,scoresPart5,scoresPart6,scoresPart7)

##write to Excel
write.xlsx(myList, file = "../results/triwiseResults/summary_GOterms_new2.xlsx")
write.xlsx(myList, "../results/triwiseResults/summary_GOterms_rank.xlsx")
write.xlsx(myList, "../results/triwiseResults/summary_GOterms_r.xlsx")
write.xlsx(myList, "../results/triwiseResults/summary_GOterms_diffexp.xlsx")
detach(package:openxlsx, unload=TRUE)

pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/General/Plots_GOterms_diffexp.pdf"), width=10)
for (Plot in 1:length(scorePlots)) {
  P<-plotPvalplot(scorePlots[[Plot]], colnames(expTable_mean))
  print(P)
}
dev.off()

plotPvalplot(scoresPart1, colnames(expTable_mean))
plotPvalplot(scoresPart2, colnames(expTable_mean))
plotPvalplot(scoresPart3, colnames(expTable_mean))
plotPvalplot(scoresPart4, colnames(expTable_mean))
plotPvalplot(scoresPart5, colnames(expTable_mean))
plotPvalplot(scoresPart6, colnames(expTable_mean))
plotPvalplot(scoresPart7, colnames(expTable_mean))

### Understand the 'angles'
scoresTest<-tail(scores,7)
scoresTest$angle<-0:6
scoresTest$qval[2]<-0.0001
plotPvalplot(scoresTest, colnames(expTable_mean))

### Plot all GO terms in specific corner
scoresOfInterest<-filter(scores, angle > 4.5, angle < 5.5)
dim(scoresOfInterest)

plotPvalplot(scoresOfInterest, colnames(expTable_mean))

scoresOfInterest$name

gsetid<-scoresOfInterest$gsetid[20]
p1<-plotDotplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F, Coi=c("")) + 
  ggtitle(scoresOfInterest$name[20]) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)


######################################################################
########## PLOT OWN GO TERMS
######################################################################

########################################
##### UPR
########################################

# ##### List 1 ====> DECIDED TO REMOVE
# UPRgenesList1<-unique(c(gsets$`GO:0035966`,gsets$`GO:0006986`,gsets$`GO:0034620`))
# UPRgenesList1<-intersect(UPRgenesList1, rownames(expTable_mean))
# 
# ## GO:0035966 = response to topologically incorrect protein
# ## GO:0006986 = response to unfolded protein
# ## GO:0034620 = cellular response to unfolded protein
# 
# p1<-plotDotplot(barycoords, Gdiffexp, UPRgenesList1, showlabels = F, Coi=c("")) + 
#   ggtitle("UPRgenesList1") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=UPRgenesList1, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)

##### List 2
UPRgenesList2<-unique(c(gsets$`GO:1990440`,gsets$`GO:0030968`,gsets$`GO:0034620`,gsets$`GO:0034976`,gsets$`GO:0006986`))
nonExpGenes <-setdiff(UPRgenesList2, rownames(expTable_mean))
UPRgenesList2<-intersect(UPRgenesList2, rownames(expTable_mean))

## GO:1990440 = positive regulation of transcription from RNA polymerase II promoter in response to endoplasmic reticulum stress
## GO:0030968 = endoplasmic reticulum unfolded protein response
## GO:0034620 = cellular response to unfolded protein
## GO:0034976 = response to endoplasmic reticulum stress
## GO:0006986 = response to unfolded protein

p1<-plotDotplot(barycoords, Gdiffexp, Goi=UPRgenesList2, showlabels = F, Coi=c("")) + 
  ggtitle("UPRgenesList2") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=UPRgenesList2, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)


### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=UPRgenesList2, showlabels = F) +
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_UPR_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=UPRgenesList2, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_UPR_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(UPRgenesList2, Gdiffexp)
nonSigGenes<-setdiff(UPRgenesList2, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "UPR_DE")
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "UPR_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "UPR_nonExp", append = T)
}

### Check how the 'angle' works
tmp<-barycoords
myGenes<-head(sigGenes,7)
tmp[myGenes,3]<-c(0,1,2,3,-3,-2,-1)
p1<-plotDotplot(tmp, Gdiffexp, Goi=myGenes[4], showlabels = F, Coi=c("")) + 
  ggtitle("UPRgenesList1") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p1)


########################################
##### UPR: Xbp1
########################################

##### Literatuur acosta-alvear
Xbp1Genes_v1<-c("Alg12","Ap3m1","Apba3","Arcn1","Arf4","Arfip2","Asb3","Bet1","Bet1l","Cd68","Cdk2ap2","Cdk5rap3","Cdk7","Copa","Copb2",
                "Copg","Creb3l2","Creld2","Dad1","Ddit3","Derl1","Dhdds","Dnajb9","Dnajb11","Dnpep","Edem1","Edem2","Eif2ak3","Elovl2",
                "Entpd7","Eny2","Ero1lb","Fam18b","Fkbp2","Gabarap","Gdap2","Ggcx","Gla","Gle1 l","Golga1","Golga3","Gys1","H13","H47",
                "Hbp1","Hdldp","Herpud1","Hist1h2ad","Hmbox1","Ift20","Impdh1","Isg20","Kbtbd3","Lman2","Lrrc59","Mbnl2","Mdh1","Mef2c",
                "Napa","Nbr1","Ncstn","Nfxl1","Oat","P4hb","Pde4dip","Pdia3","Pdia6","Pecr","Pgm3","Pla2g6","Polr3k","Ppib","Ppp1r10",
                "Prdx6","Preb","Rpn1","Sar1a","Sec13l1","Sec23b","Sec24d","Sec63","Selm","Slc2a6","Slc30a5","Slc30a7","Spcs2","Srp19",
                "Srp54b","Srpr","Srprb","Ssr2","Ssr3","Stch","Stt3a","Stx5a","Stx18","Tbc1d15","Tipf5","Tme3","Tmem39a","Tmem167b",
                "Tnfaip1","Txndc4","Txndc11","Ube1dc1","Ufd1l","Ufm1","Yipf2","Zfand2a")
Xbp1Genes_v1<-c(Xbp1Genes_v1,"Txndc5", "Dnajc10", "Sec61a1") #Adjustment Sophie summer 2019
setdiff(Xbp1Genes_v1, rownames(expTable_mean))
setdiff(Xbp1Genes_v1, rownames(countData))
Xbp1Genes_v1[Xbp1Genes_v1=="Copg"]<-"Copg1"
Xbp1Genes_v1[Xbp1Genes_v1=="Fam18b"]<-"Tvp23b"
Xbp1Genes_v1[Xbp1Genes_v1=="Gle1 l"]<-"Gle1"
Xbp1Genes_v1[Xbp1Genes_v1=="H47"]<-"Selenos"
Xbp1Genes_v1[Xbp1Genes_v1=="Sec13l1"]<-"Sec13"
Xbp1Genes_v1[Xbp1Genes_v1=="Stch"]<-"Hspa13"
Xbp1Genes_v1[Xbp1Genes_v1=="Txndc4"]<-"Erp44"
Xbp1Genes_v1[Xbp1Genes_v1=="Ube1dc1"]<-"Uba5"
nonExpGenes <-setdiff(Xbp1Genes_v1, rownames(expTable_mean))
Xbp1Genes_v1<-intersect(Xbp1Genes_v1, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v1, showlabels = F, Coi=c("")) + 
  ggtitle("Xbp1 literatuur") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v1, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v1, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_UPRXbp1_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v1, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_UPRXbp1_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(Xbp1Genes_v1, Gdiffexp)
nonSigGenes<-setdiff(Xbp1Genes_v1, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Xbp1_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Xbp1_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Xbp1_nonExp", append = T)
}

# ##### GO term IRE1 GO:0036489, downloaded from Amigo (Human) ====> DECIDED TO REMOVE
# Xbp1Genes_v2<-c("ERN1","XBP1","PTPN1","HSPA5","LMNA","BAX","DAB2IP","HSPA1A","DCTN1","WFS1","DDX11","BCL2L11","BAK1","TPP1","SHC1",
#                 "ACADVL","COPS5","GSK3A","SYVN1","SEC31A","ADD1","SULT1A3","CUL7","WIPI1","PDIA6","PLA2G4B","TMBIM6","ZBTB17","DNAJC3",
#                 "AGR2","GOSR2","EXTL3","TLN1","GFPT1","CXXC1","DNAJB11","ASNA1","HDGF","EXTL2","TSPYL2","HYOU1","ATP6V0D1","BFAR",
#                 "SEC63","SEC61B","PPP2R5B","TMEM33","CTDSP2","EDEM1","SEC61A1","MYDGF","ARFGAP","SSR1","SERP1","PREB","PDIA5","SRPRA",
#                 "YIF1A","FKBP14","SEC62","DNAJB9","EXTL1","SRPRB","KDELR3","SEC61G","KLHDC3","SEC61A2","TATDN2")
# Xbp1Genes_v2<-firstup(tolower(Xbp1Genes_v2))
# setdiff(Xbp1Genes_v2, rownames(expTable_mean))
# setdiff(Xbp1Genes_v2, rownames(countData))
# Xbp1Genes_v2<-intersect(Xbp1Genes_v2, rownames(expTable_mean))
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v2, showlabels = F, Coi=c("")) + 
#   ggtitle("Xbp1 GO term") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Xbp1Genes_v2, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)
# 
# 
# ##### Xbp1 merge ====> DECIDED TO REMOVE
# Xbp1Genes<-unique(c(Xbp1Genes_v1, Xbp1Genes_v2))
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=Xbp1Genes, showlabels = F, Coi=c("")) + 
#   ggtitle("Xbp1 merge") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Xbp1Genes, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)


########################################
##### UPR: Isr
########################################

# ##### Common targets ====> DECIDED TO REMOVE
IsrCommon<-c("Atf3","Bcat2","Bdnf","Cebpg","Clcn3","Cxadr","Ddit3","Eif4g2","B4galnt2","Gng5","Hspa5","Has2","Hfe","Cd74","Incenp",
             "Itgb5","Itgb7","Itih1","Kif13b","Kifc1","Arhgef2","Lcn2","Lgals3","Amacr","Aff1","Mthfd2","Mtm1","Ppp1r15a","Nfe2l1",
             "Nme2","Ddr2","Sqstm1","Pax8","Per2","Pik3c2g","Prkg2","Ptbp1","Ptn","Pvrl2","Reln","Dpf2","Rpl7","Sars","Slc20a1",
             "Slc7a5","Snai2","Stc2","Tagln2","Tgif1","Tnks","Rpl13a","Ube2l3","Ube2g2","Uqcrq","Vars","Vldlr","Wars","Ywhag",
             "Ets1","Pacsin2","Ptpn21","Slc7a11","Deb1","Dysf","Zfp238","Ero1l","D4Bwg0951e","Pnrc2","Gnpnat1","Gtpbp2","Htra1",
             "Tspan5","Eif3c","Aldh18a1","Cnpy2","Ylpm1","Ick","Nfu1","Sertad2","Trem3","Nrip2","Fads3","Fam129a","B230120H23Rik",
             "Krtcap2","Ghitm","Nosip","2610528E23Rik","Ube2w","Wwp2","Edem3","Lsm14a","Naaa","Pfdn1","Eif2s2","Dusp6","Fibin",
             "Ift172","Mrpl24","Slc25a39","Otub2","Leprotl1","Fam96a","Lmbrd1","Tab2","Flnc","1810031K17Rik","Degs2","Nars","Hdac8",
             "Steap1","2810408A11Rik","Taf15","Btf3l4","Fam175a","Angptl6","Slc41a3","Coasy","Cyb5r1","Tpx2","Kbtbd5","Tspyl4",
             "2810006K23Rik","Golim4","Xpot","Rd3","Lonp1","1300001I01Rik","Stk40","1700016K19Rik","4931414P19Rik","Samd4","Ascc2",
             "Zc3h18","Snrnp35","Tbc1d9b","Zc3hav1","Rnf114","Il23a","Dnaja3","Rbm9","Ugt1a6a","Hist2h4","Usp6nl","Akna","Emilin1",
             "6330503K22Rik","Iars","Mycbp2","Lars","Macrod1","Yars","Atf5","Eprs","Ankrd1","Slco1a5","Gpt2","Rcc2","Ppp1r15b","Dut",
             "Srsf1","Acaa1a","Prdm15","Slc19a2","Pop5","Rbm39","Acot2","Zbtb7c","2700007P21Rik","Zfp598","Fam83f","A630007B06Rik",
             "Aldh1l2","Zfc3h1","Mars","Tmem11","Kdm6b","Eif5","Lhfpl2","March6","Soat2","Cyp2ab1","Ubr2","Erlin1","Atf6","Camsap1",
             "Phyhd1","Trib3","Rsbn1","Paqr3","Lrrk1","Cnot1","Slc38a7","Aars","Snx33","Tbpl1","Cdh24","Dleu7","Abcc4","A630033E08Rik",
             "Slc1a7","Gpatch3","Zbtb38","Pkmyt1","Optc","Tmtc2","Gm867","Jhdm1d","Agap1","Gars","Ypel5","Cdsn","Flrt1","Tmem189",
             "Ndufa4l2","Fam159a","Gpx4","Defb43")
setdiff(IsrCommon, rownames(expTable_mean))
setdiff(IsrCommon, rownames(countData))
IsrCommon[IsrCommon=="Zfp238"]<-"Zbtb18"
IsrCommon[IsrCommon=="D4Bwg0951e"]<-"Lurap1l"
IsrCommon[IsrCommon=="B230120H23Rik"]<-"Map3k20"
IsrCommon[IsrCommon=="2610528E23Rik"]<-"Cmss1"
IsrCommon[IsrCommon=="1810031K17Rik"]<-"Cnppd1"
IsrCommon[IsrCommon=="Kbtbd5"]<-"Klhl40"
IsrCommon[IsrCommon=="1300001I01Rik"]<-"Cluh"
IsrCommon[IsrCommon=="Rbm9"]<-"Rbfox2"
IsrCommon[IsrCommon=="6330503K22Rik"]<-"Ccp110"
IsrCommon[IsrCommon=="2700007P21Rik"]<-"Arl14ep"
IsrCommon[IsrCommon=="A630007B06Rik"]<-"Ccdc186"
IsrCommon[IsrCommon=="A630033E08Rik"]<-"Zfp945"
IsrCommon[IsrCommon=="Jhdm1d"]<-"Kdm7a"
nonExpGenes <-setdiff(IsrCommon, rownames(expTable_mean))
IsrCommon<-intersect(IsrCommon, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=IsrCommon, showlabels = F, Coi=c("")) +
  ggtitle("Isr common") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrCommon, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=IsrCommon, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_IsrCommon.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrCommon, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_IsrCommon.png",dpi=300)


##### CHOP targets ====> DECIDED TO REMOVE
# IsrChop<-c("Abl1","Akr1b7","Car6","Cav1","Entpd5","Csrp1","Emp2","Eps8","Erh","Gpc4","Hist2h3c1","Id1","Insr","Itih2","Jak2",
#            "Mycl1","Tm4sf1","Mbd2","Mbd4","Mcm2","Mds1","Me1","Msmb","Mt2","Mttp","Ngf","Nr4a2","Pafah1b2","Psmc2","Ptpn2","C80913",
#            "Snrk","Star","Trim24","Trhr","Ywhaz","Nt5e","Pibf1","Sdcbp","Tacc2","Arl6ip5","Llph","Zfp869","Wls","Fam135a","Rpa3",
#            "4930547C10Rik","Ms4a6d","Krtap3-1","2510012J08Rik","Wars2","Pwwp2a","9130011E15Rik","Ccdc19","Dus4l","Actr3","Ttc27",
#            "Gzf1","4930486G11Rik","Ndfip2","1700081L11Rik","Lass2","Recql4","Dach2","Klra18","Nol12","Wdr92","Kazald1","Camk2d",
#            "Fam64a","Lims1","Palmd","Ehbp1l1","Il17rd","Zfp184","Osgin2","Vmn1r183","Tmem60","D4Ertd22e","Arhgef11","Sipa1l1",
#            "Chdh","Reep2","Ralgapb","Hnrnpa3","Uba6","Agpat9","Foxred1","Ostn","Dusp5","Tnfsf18","Rpusd2","Nrcam","Islr2","Bend6",
#            "Morc3","Dgkh","Parl","9630014M24Rik","Defb45","Znhit3","Kncn","Vmn1r223")
# setdiff(IsrChop, rownames(expTable_mean))
# setdiff(IsrChop, rownames(countData))
# IsrChop[IsrChop=="Mycl1"]<-"Mycl"
# IsrChop[IsrChop=="Mds1"]<-"Mecom"
# IsrChop[IsrChop=="C80913"]<-"Uri1"
# IsrChop[IsrChop=="4930547C10Rik"]<-"Toporsl"
# IsrChop[IsrChop=="2510012J08Rik"]<-"Cactin"
# IsrChop[IsrChop=="Ccdc19"]<-"Cfap45"
# IsrChop[IsrChop=="4930486G11Rik"]<-"Mei4"
# IsrChop[IsrChop=="1700081L11Rik"]<-"Kansl1"
# IsrChop[IsrChop=="Lass2"]<-"Cers2"
# IsrChop[IsrChop=="D4Ertd22e"]<-"Szrd1"
# IsrChop<-intersect(IsrChop, rownames(expTable_mean))
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=IsrChop, showlabels = F, Coi=c("")) +
#   ggtitle("Isr Chop") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrChop, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)
# 
# ### Export plots
# p<-plotDotplot(barycoords, Gdiffexp, Goi=IsrChop, showlabels = F) + 
#   theme(legend.position = "none")
# ggsave(p,file="../results/triwiseResults/plots/triwisePlot_IsrCHOP.png",dpi=300)
# 
# p<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrChop, showlabels = F)
# ggsave(p,file="../results/triwiseResults/plots/rosePlot_IsrCHOP.png",dpi=300)

##### ATF4 targets
IsrAtf<-c("Akap2","Apbb2","Apobec1","Slc7a1","Bcat1","Casp4","Cast","Runx3","Col18a1","Cp","Bcar1","Dpysl2","Cr1l","Ddc","Eif4ebp1",
          "F3","Fcer1g","Tsc22d3","Slc6a9","Got1","Gpam","Grb10","Gsto1","H2-Q7","Hmox1","Ier2","Il1r1","Il2","Inpp5b","Mdfic",
          "Klf4","Serpina3c","Lmo4","Ltbp3","Maz","Anapc1","Slc3a2","Mras","Mtap4","Mthfr","Myom1","Ncl","Ndrg1","Nf1","Nfkb1",
          "Nos2","Osmr","Pde1a","Per3","Abcb1a","Pld1","Ptpn12","Ptprs","Rad21","Rai1","Ran","Rbm4","Rps6ka2","Ruvbl2","Ccl2","Ccl7",
          "Serpinf1","Slc1a5","Slit1","Serpinb9","Srpk1","Eif1","Surf6","Synj2","Tnfaip2","Uba3","Ubp1","Wbp4","Wisp1","Zbtb7b",
          "Hax1","Vnn3","Asns","Cars","Cabp1","Clic4","Tor3a","H2afz","Usp2","Htatip2","Tslp","Calcrl","Slc1a4","Snx12","Lgals8",
          "Atl2","Nupr1","Spen","Fgf21","Dnajb12","Angptl4","Mmp19","Slc35b4","Mettl9","Bhlhe22","Bcmo1","Herpud1","Gpr85","Mrpl54",
          "1110018J18Rik","Klf15","Atp6v1g1","Pla2g12a","Dhrs7","Crls1","Tmigd1","Ormdl3","Sltm","4921524J17Rik","Trim35","Acad8",
          "Uqcrc2","Ttc23","Fam119a","2900073G15Rik","3110049J23Rik","Gcfc1","Anp32b","Samd8","Slc38a2","Jagn1","Rab39b","Mid1ip1",
          "B230217C12Rik","Rnaset2b","Mrpl13","Chac1","Ubr4","Bola1","Gtpbp4","Ipmk","Ttc9c","Tbce","Slc25a33","Nipbl","Arid5b",
          "Wdr27","Cul2","Slc16a14","Phf10","Fbxl2","Usp36","Spats2","Arpc5l","Pck2","Scpep1","Ddit4","Rprd2","Mff","Fbxo31","Casc5",
          "Ip6k2","Klhdc10","Rhbdd1","Ankrd11","Prss36","Pgpep1l","Zfp623","Crispld2","Sh3bp5l","Jdp2","Pde4dip","Acox2","Echs1",
          "Dtnbp1","Wwtr1","Ifi44","Cd276","Traf3ip2","Rhoq","Dhx57","Psat1","Acaca","Cth","Shmt2","Mat2b","Fam175b","Ank2","Tars",
          "Spred1","Stard5","Phf21a","Tet3","Daam1","Pycr1","Leprel1","Cln5","Ifitm6","Hhipl1","Kdm5a","Sema6d","Fam26f","Zfp365",
          "Plekhh3","Tnrc6c","Khnyn","Adm2","2310008H04Rik","Lsg1","Crybg3","Bnip1","Uhrf1bp1","Fbxo11","Etf1","Taf5","Ppp2r5a",
          "Gm129","Sesn2","Pdap1","B9d2","Atxn2l","BC019943","Inpp4b","Vac14","Phgdh","Usp11","B230315N10Rik","C130039O16Rik",
          "Tnpo1","Tmem74","Kcnt2","Manea","BC055111","Cct8l1","Prss35","Vps54","Slc25a28","Alkbh5","Zfp608","Mthfd1l","Flad1",
          "Fndc7","Tmem154","5930403L14Rik","4933421E11Rik","Cxcl3","Fam188b","Trim66","Slc9a9","Tlcd2","Zyg11b","Dnaic2","Zfp708",
          "BB014433","Gm5627","Gm5820","Gsk3a","Nrbf2","9130023H24Rik")
setdiff(IsrAtf, rownames(expTable_mean))
setdiff(IsrAtf, rownames(countData))
IsrAtf[IsrAtf=="Mtap4"]<-"Map4"
IsrAtf[IsrAtf=="Bcmo1"]<-"Bco1"
IsrAtf[IsrAtf=="1110018J18Rik"]<-"Aaed1"
IsrAtf[IsrAtf=="Fam119a"]<-"Mettl21a"
IsrAtf[IsrAtf=="2900073G15Rik"]<-"Myl12a"
IsrAtf[IsrAtf=="3110049J23Rik"]<-"Pbld2"
IsrAtf[IsrAtf=="Gcfc1"]<-"Paxbp1"
IsrAtf[IsrAtf=="Leprel1"]<-"P3h2"
IsrAtf[IsrAtf=="2310008H04Rik"]<-"Spidr"
IsrAtf[IsrAtf=="Gm129"]<-"Ciart"
IsrAtf[IsrAtf=="BC019943"]<-"Tti2"
IsrAtf[IsrAtf=="B230315N10Rik"]<-"Zfp938"
IsrAtf[IsrAtf=="C130039O16Rik"]<-"Elmsan1"
IsrAtf[IsrAtf=="4933421E11Rik"]<-"Lrif1"
nonExpGenes <-c(nonExpGenes, setdiff(IsrAtf, rownames(expTable_mean))) #Add to list of nonEXPgenes from IsrCommon
IsrAtf<-intersect(IsrAtf, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=IsrAtf, showlabels = F, Coi=c("")) + 
  ggtitle("Isr Atf4") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrAtf, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=IsrAtf, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_IsrAtf_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrAtf, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_IsrAtf_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(IsrAtf, Gdiffexp)
nonSigGenes<-setdiff(IsrAtf, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
# write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "IsrAtf4_DE", append = T)
# write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "IsrAtf4_nonDE", append = T)
# if (length(nonExpGenes) > 1) {
#   write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "IsrAtf4_nonExp", append = T)
# }

# ##### Common + ATF4 targets ====> DECIDED TO REMOVE
Isr<-unique(c(IsrCommon, IsrAtf))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Isr, showlabels = F, Coi=c("")) +
  ggtitle("Isr") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Isr, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Isr, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_IsrCommon_and_Atf.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Isr, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_IsrCommon_and_Atf.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(Isr, Gdiffexp)
nonSigGenes<-setdiff(Isr, Gdiffexp)
nonExpGenes <-unique(nonExpGenes)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/General/summary_ownGOlists_final_paper_oct2019.xls", "Isr_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/General/summary_ownGOlists_final_paper_oct2019.xls", "Isr_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/General/summary_ownGOlists_final_paper_oct2019.xls", "Isr_nonExp", append = T)
}
# 
# ##### GO terms ====> DECIDED TO REMOVE
# IsrGOgenes<-unique(c(gsets$`GO:0006520`,gsets$`GO:00070059`))
# setdiff(IsrGOgenes, rownames(expTable_mean))
# setdiff(IsrGOgenes, rownames(countData))
# IsrGOgenes<-intersect(IsrGOgenes, rownames(expTable_mean))
# 
# ## GO:0006520 = cellular amino acid metabolism
# ## GO:00070059 = intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=IsrGOgenes, showlabels = F, Coi=c("")) + 
#   ggtitle("Isr GO terms") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=IsrGOgenes, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)


########################################
##### UPR: Atf6
########################################

atf6<-c("Bip","Grp94","Pdia4","Ero1l","Crt","Herpud1","Os9","Sel1l","Erdj3","Erp57","Pdia6","Vcp","Uggt1","Sec11c","P58","Hyou1",
        "Pdia10","Edem1","Derl2","Sec13")
setdiff(atf6, rownames(expTable_mean))
setdiff(atf6, rownames(countData))
atf6[atf6=="Bip"]<-"Hspa5"
atf6[atf6=="Grp94"]<-"Hsp90b1"
atf6[atf6=="Crt"]<-"Slc6a8"
atf6[atf6=="Erdj3"]<-"Dnajb11"
atf6[atf6=="P58"]<-"Dnajc3"
atf6[atf6=="Erp57"]<-"Pdia3"
nonExpGenes <-setdiff(atf6, rownames(expTable_mean))
atf6<-intersect(atf6, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=atf6, showlabels = F, Coi=c("")) + 
  ggtitle("Atf6") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=atf6, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=atf6, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_Atf6_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=atf6, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_Atf6_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(atf6, Gdiffexp)
nonSigGenes<-setdiff(atf6, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Atf6_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Atf6_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Atf6_nonExp", append = T)
}

########################################
##### RIDD
########################################

##### Literatuur Han et al, 2009
riddHan<-c("DNAJC3","TXNDC5","NPC1","GUSB","DGAT1","SLC36A4","M6PR","GNPTG","CSTF2T","TMEM131","MGAT4B","PLXNA1","NXPH4","PXDN","APLP1",
           "TMEM129","FRAS1","CDS2","CTSB","CBARA1","LEPRE1","TPP1","FAM8A1","ITM2C","ASAH1","CMTM8","CPE","CPD","PLD3","ADAM9","NDUFS8",
           "ZDHHC24","SUMF2","APP","PTPRG","GTF3C1","ATP2B4","NPC2","B2M","SEMA6A","LAMC1","PAQR3","SERINC3","ERGIC2","KTN1","NEO1","TAP1",
           "ZDHHC9","TSG101","TSPAN33","TMEM64","SCPEP1","RYK","GSTK1","SLC39A14","LRRC8C","TAPBP","CNNM3","BAT5","HEXB","FURIN","KAL1",
           "PLXNB2","C9orf5","SH3GL2","ENPP2","RCN1","TWSG1","GALNT1","PRAF2","GLB1","JAGN1","NOMO2","ATP1A1","TM9SF2","TSPAN3","SORT1",
           "SPPL2A","ADPGK","REEP3","GPAA1","DNAJC1","RNF167","ITGB1","TXNDC12","GGH","TUSC3","TEX264","RCN2","PCDH7","RNF130","NIPA2",
           "PRKCSH","LAMP1","ADAM15","CTSC","KDELR2","SLC35B1","EPHX1","KDELC1","ZDHHC16","TMEM48","IGFBP2","PPP2R1A","CD151","LAMA4",
           "SCAMP2","PRSS23","COMT","YIPF3","ABCC4","TOMM20","SLC39A1","CMTM6","ATP9A","ACSL3","CDIPT","RQCD1","PRDX4","NPTX2","SQLE",
           "HSD17B12","OCIAD1","LDLRAD3","LGR4","CALR","ATP2A2","CKAP4","TMED9","APOA1BP","SIAE","LRP11","ATP2B1","SLIT2","AMFR","ASPH",
           "EIF2AK1","CST3","TMED2","LMAN2","TMEM38B","REEP2","SULF2","XK","GPC4","DSG2","TMED5","CHPT1","MXRA7","BACE2","CASC4","TUFM",
           "RPN2","NOC2L","PDE3B","TMEM50A","AGPAT2","HMGCR","ANTXR1","ICMT","PTTG1IP","DSC2","PTPRK","SGCE","ITPR1","PVRL2","TGFBR3",
           "CD99L2","LAMB1","SEC63","FAM3A","LDOC1","STOM","C6orf192","CD9","TMEM132A","TMTC1","SLC25A1","PTPRF","CMTM4","TMCO3","PRKD1",
           "PTPLA","RPAP1","SLC38A6","DDR1","SLC39A10","PPAPDC2","SLC6A15","CBS","NQO2","HTATIP2","B4GALT1","COLEC12","IFI30","LMAN1",
           "PPIB","PPIC","FADS3","DONSON")
riddHan<-firstup(tolower(riddHan))
setdiff(riddHan, rownames(expTable_mean))
setdiff(riddHan, rownames(countData))
riddHan[riddHan=="Cbara1"]<-"Micu1"
riddHan[riddHan=="Lepre1"]<-"P3h1"
riddHan[riddHan=="Bat5"]<-"Abhd16a"
riddHan[riddHan=="Kal1"]<-"Slamf6"
riddHan[riddHan=="Tmem48"]<-"Ndc1"
riddHan[riddHan=="Ptpla"]<-"Hacd1"
riddHan[riddHan=="C9orf5"]<-"Tmem245"
riddHan[riddHan=="C6orf192"]<-"Slc18b1"
nonExpGenes <-setdiff(riddHan, rownames(expTable_mean))
nonExpGenes2 <-setdiff(riddHan, rownames(expTable_mean)) #For merged list
riddHan<-intersect(riddHan, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=riddHan, showlabels = F, Coi=c("")) + 
  ggtitle("Ridd Han et al, 2009") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=riddHan, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=riddHan, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_riddHan_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=riddHan, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_riddHan_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(riddHan, Gdiffexp)
nonSigGenes<-setdiff(riddHan, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "RiddHan_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "RiddHan_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "RiddHan_nonExp", append = T)
}

# ##### Literatuur Sarah Gerlo 
riddGerlo<-c("ANGPTL3","BIP","BLOC1","CD59","CES1","COL6A1","CYP1A2","CYP2E1","ERGIC3","GALNT2","GEMIN5","GPC3","GYLTL1B",
            "HGSNAT","INS","IRE1a","ITGB2","MIR17","MKRN2","PDGFRB","PDIA4","PDK2","PEPD","PER1","PPP2R1A","PMP22","PR-4","PRKCD",
            "RTN4","RUVBL1","SCARA3","SPARC","SUMO","28S","TAPBP","ms","YWHAQ")
riddGerlo<-firstup(tolower(riddGerlo))
setdiff(riddGerlo, rownames(expTable_mean))
setdiff(riddGerlo, rownames(countData))
riddGerlo[riddGerlo=="Bip"]<-"Hspa5"
riddGerlo[riddGerlo=="Cd59"]<-"Cd59a"
riddGerlo[riddGerlo=="Ces1"]<-"Ces1g"
riddGerlo[riddGerlo=="Ire1a"]<-"Ern1"
riddGerlo[riddGerlo=="28s"]<-"Rn28s1"
riddGerlo[riddGerlo=="Ms"]<-"Mtr"
riddGerlo[riddGerlo=="Bloc1"]<-"Bloc1s1" #From Sophie. Alias not found on NCBI.
nonExpGenes <-setdiff(riddGerlo, rownames(expTable_mean))
nonExpGenes2 <-c(nonExpGenes2, setdiff(riddGerlo, rownames(expTable_mean))) #For merged list
riddGerlo<-intersect(riddGerlo, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=riddGerlo, showlabels = F, Coi=c("")) +
 ggtitle("Ridd Sarah Gerlo") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=riddGerlo, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=riddGerlo, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_riddGerlo_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=riddGerlo, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_riddGerlo_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(riddGerlo, Gdiffexp)
nonSigGenes<-setdiff(riddGerlo, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "riddGerlo_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "riddGerlo_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "riddGerlo_nonExp", append = T)
}
 
##### Ridd merge 
riddMerge<-unique(c(riddHan, riddGerlo))
setdiff(riddMerge, rownames(expTable_mean))
setdiff(riddMerge, rownames(countData))
nonExpGenes2 <-unique(nonExpGenes2) #Merged list
riddMerge<-intersect(riddMerge, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=riddMerge, showlabels = F, Coi=c("")) +
 ggtitle("Ridd merge") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=riddMerge, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=riddMerge, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_riddMerge_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=riddMerge, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_riddMerge_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(riddMerge, Gdiffexp)
nonSigGenes<-setdiff(riddMerge, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "RiddMerge_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "RiddMerge_nonDE", append = T)
if (length(nonExpGenes2) > 1) {
  write.xlsx(nonExpGenes2, "../results/heatmaps/Final_paper_Simon/General/summary_ownGOlists_final_paper_oct2019.xls", "RiddMerge_nonExp", append = T)
}


########################################
##### Lipid metabolism
########################################

##### All lipid
lipidGenes<-unique(c(gsets$`GO:1903727`,gsets$`GO:1903725`,gsets$`GO:0045834`,gsets$`GO:0046890`,gsets$`GO:0019216`,gsets$`GO:0006629`,
                     gsets$`GO:0090181`,gsets$`GO:0050810`,gsets$`GO:0019218`))
setdiff(lipidGenes, rownames(expTable_mean))
setdiff(lipidGenes, rownames(countData))
nonExpGenes <-setdiff(lipidGenes, rownames(expTable_mean))
lipidGenes<-intersect(lipidGenes, rownames(expTable_mean))

## GO:1903727 = positive regulation of phospholipid metabolic process
## GO:1903725 = regulation of phospholipid metabolic process
## GO:0045834 = positive regulation of lipid metabolic process
## GO:0046890 = regulation of lipid biosynthetic process
## GO:0019216 = regulation of lipid metabolic process
## GO:0006629 = lipid metabolic process
## GO:0090181 = regulation of cholesterol metabolic process
## GO:0050810 = regulation of steroid biosynthetic process
## GO:0019218 = regulation of steroid metabolic process

p1<-plotDotplot(barycoords, Gdiffexp, Goi=lipidGenes, showlabels = F, Coi=c("")) + 
  ggtitle("lipidGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=lipidGenes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=lipidGenes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_lipid_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=lipidGenes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_lipid_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(lipidGenes, Gdiffexp)
nonSigGenes<-setdiff(lipidGenes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Lipid_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Lipid_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Lipid_nonExp", append = T)
}

########################################
##### Cholesterol: metabolism
########################################
humanGenes<-c("ACAT2","APOB","FDPS","FDPS","HMGCS2","PRKAA1","PRKAA1","ACACA","HMGCR","HMGCR","HMGCR","ERLIN2","PLPP6","CES1","IDI1",
              "IDI1","NFYC","EBP","EBP","EBP","SP1","GPAM","POR","NFYB","SREBF2","SCAP","SCAP","LSS","LSS","LSS","MVK","MVK","MVK",
              "INSIG1","PMVK","PMVK","PMVK","IDI2","DHCR24","DHCR24","DHCR24","DHCR24","MSMO1","SREBF1","SREBF1","TM7SF2","TM7SF2",
              "ACLY","NSDHL","MVD","MVD","MVD","MVD","G6PD","APOA4","ABCG1","LBR","MBTPS1","FGF1","PRKAA2","CNBP","SC5D","SC5D","SC5D",
              "SC5D","NPC1L1","CYB5R3","ACAA2","APOE","APOA1","HMGCS1","HMGCS1","ARV1","DHCR7","DHCR7","DHCR7","DHCR7","CYP51A1",
              "CYP51A1","CYP51A1","FASN","SOD1","NFYA","SCD","ACACB","INSIG2","APOA5","CFTR","HSD17B7","SEC14L2","ERLIN1","MBTPS2",
              "RAN","GGPS1","GGPS1","KPNB1","ELOVL6","FDFT1","FDFT1","FDFT1","SQLE","SQLE")
mouseGenes<-c("Erlin1","Hmgcs2","Pmvk","Pmvk","Fgf1","Mvk","Mvk","Scap","Hmgcs1","Cftr","Fdft1","Fdft1","Cyp51","Cyp51","Cyp51","Apoe",
              "Apoe","Apoe","Apob","Apoa4","Apoa1","Nsdhl","Scp2","Ebp","Ebp","Abcg1","Prkaa2","Lss","Lss","Npc1l1","Insig2","Hmgcr",
              "Hmgcr","Hmgcr","Hmgcr","Hmgcr","Erlin2","Cyb5r3","Fdps","Fdps","Prkaa1","Insig1","Sc5d","Sc5d","Hsd17b7","Idi2","Mvd",
              "Dhcr7","Dhcr7","Dhcr7","Dhcr7","G6pdx","Pex2","Sod1","Sec14l2","Apoa5","Por","Ces1d","Tm7sf2","Dhcr24","Dhcr24","Dhcr24",
              "Dhcr24","Srebf1","Srebf1")
cholesterolGenes<-unique(firstup(tolower(c(humanGenes,mouseGenes))))
cholesterolGenes<-c(cholesterolGenes,"Nus1","Lss","Tm7sf2","Sc4mol","Nsdhl","Dhcr24","Sc5d","Dhcr7","Ebp") #Adjustment Sophie
setdiff(cholesterolGenes, rownames(expTable_mean))
setdiff(cholesterolGenes, rownames(countData))
cholesterolGenes[cholesterolGenes=="Plpp6"]<-"Ppapdc2"
cholesterolGenes[cholesterolGenes=="Ces1"]<-"Ces1g"
cholesterolGenes[cholesterolGenes=="G6pd"]<-"G6pdx"
cholesterolGenes[cholesterolGenes=="Scd"]<-"Scd1"
cholesterolGenes[cholesterolGenes=="Cyp51a1"]<-"Cyp51"
cholesterolGenes<-unique(cholesterolGenes)
nonExpGenes <-setdiff(cholesterolGenes, rownames(expTable_mean))
cholesterolGenes<-intersect(cholesterolGenes, rownames(expTable_mean))


p1<-plotDotplot(barycoords, Gdiffexp, Goi=cholesterolGenes, showlabels = F, Coi=c("")) + 
  ggtitle("cholesterolGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=cholesterolGenes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=cholesterolGenes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_cholesterolMetabolism_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=cholesterolGenes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_cholesterolMetabolism_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(cholesterolGenes, Gdiffexp)
nonSigGenes<-setdiff(cholesterolGenes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_met_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_met_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_met_nonExp", append = T)
}

########################################
##### Cholesterol: efflux
########################################
effluxGenes<-unique(c(gsets$`GO:0033344`,gsets$`GO:0010874`,gsets$`GO:0010875`))
setdiff(effluxGenes, rownames(expTable_mean))
setdiff(effluxGenes, rownames(countData))
nonExpGenes <-setdiff(effluxGenes, rownames(expTable_mean))
effluxGenes<-intersect(effluxGenes, rownames(expTable_mean))


p1<-plotDotplot(barycoords, Gdiffexp, Goi=effluxGenes, showlabels = F, Coi=c("")) + 
  ggtitle("Cholesterol effluxGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=effluxGenes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=effluxGenes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_cholesterolEfflux_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=effluxGenes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_cholesterolEfflux_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(effluxGenes, Gdiffexp)
nonSigGenes<-setdiff(effluxGenes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_efflux_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_efflux_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_efflux_nonExp", append = T)
}

########################################
##### Cholesterol: biosynthese
########################################
biosyntheseGenes<-unique(c(gsets$`GO:0006695`,gsets$`GO:0045540`,gsets$`GO:0045542`))
setdiff(biosyntheseGenes, rownames(expTable_mean))
setdiff(biosyntheseGenes, rownames(countData))
nonExpGenes <-setdiff(biosyntheseGenes, rownames(expTable_mean))
biosyntheseGenes<-intersect(biosyntheseGenes, rownames(expTable_mean))


p1<-plotDotplot(barycoords, Gdiffexp, Goi=biosyntheseGenes, showlabels = F, Coi=c("")) + 
  ggtitle("Cholesterol biosyntheseGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=biosyntheseGenes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=biosyntheseGenes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_cholesterolBiosynthesis_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=biosyntheseGenes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_cholesterolBiosynthesis_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(biosyntheseGenes, Gdiffexp)
nonSigGenes<-setdiff(biosyntheseGenes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_biosynthese_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_biosynthese_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_biosynthese_nonExp", append = T)
}

########################################
##### Cholesterol: all
########################################
cholesterolGenesAll<-unique(c(cholesterolGenes, effluxGenes, biosyntheseGenes))
setdiff(cholesterolGenesAll, rownames(expTable_mean))
setdiff(cholesterolGenesAll, rownames(countData))


p1<-plotDotplot(barycoords, Gdiffexp, Goi=cholesterolGenesAll, showlabels = F, Coi=c("")) + 
  ggtitle("Cholesterol all") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=cholesterolGenesAll, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=cholesterolGenesAll, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_cholesterolAll_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=cholesterolGenesAll, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_cholesterolAll_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(cholesterolGenesAll, Gdiffexp)
nonSigGenes<-setdiff(cholesterolGenesAll, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_all_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_all_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Cholesterol_all_nonExp", append = T)
}

########################################
##### Cholesterol: AMIGO genes (Clint)
########################################

#Read in Amigo results from excel file
Amigo_chol_genes <- read_xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/documentation/Amigo_chol_genes_Sophie_04-07-19.xlsx", col_names = T)

#Create string and counter
Amigo_chol_genes_known <- ""
c = 1

#Loop over Amigo results and filter out those with known genes
for (i in 1:length(Amigo_chol_genes$Cholesterol_genes_AMIGO)) {
  matched <- str_detect(Amigo_chol_genes$Cholesterol_genes_AMIGO[i], "UniProtKB") #Lines with UniprotKB have known genes
  
  if (matched == T) {
    Long_name <- str_split(Amigo_chol_genes$Cholesterol_genes_AMIGO[i], "        ") #Split on spaces
    Gene <- Long_name[[1]][2] #2nd part holds the gene name
    Amigo_chol_genes_known[c] <- Gene
    c = c + 1
  }
}

Amigo_chol_genes_known<-unique(firstup(tolower(Amigo_chol_genes_known)))
setdiff(Amigo_chol_genes_known, rownames(expTable_mean))
setdiff(Amigo_chol_genes_known, rownames(countData))
#Check for aliases on ncbi for genes not in countdata table
Amigo_chol_genes_known[Amigo_chol_genes_known=="Akr1c1"]<-"Akr1c6"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Disp3"]<-"Ptchd2"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Relch"]<-"2310035C23Rik"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Pip4p1"]<-"Tmem55b"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Plpp6"]<-"Ppapdc2"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Cd24"]<-"Cd24a"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Apol2"]<-"Apol8"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Agtr1"]<-"Agtr1a"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Scd"]<-"Scd1"
Amigo_chol_genes_known[Amigo_chol_genes_known=="G6pd"]<-"G6pdx"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Ces1"]<-"Ces1g"
Amigo_chol_genes_known[Amigo_chol_genes_known=="Cyp51a1"]<-"Cyp51"
Amigo_chol_genes_known<-unique(Amigo_chol_genes_known)
nonExpGenes <-setdiff(Amigo_chol_genes_known, rownames(expTable_mean))
Amigo_chol_genes_known<-intersect(Amigo_chol_genes_known, rownames(expTable_mean))
Amigo_chol_genes_known<-setdiff(Amigo_chol_genes_known, c("Ccl3", "Nfkbia","Xbp1")) #Adjustment Sophie summer 2019 to avoid confusion and remove genes which shouldn't be in the GO term.
Amigo_chol_genes_known<-c(Amigo_chol_genes_known, "Apol10b","Apol7c") #Adjustment Sophie October 2019 (V2)
Amigo_chol_genes_known<-c(Amigo_chol_genes_known, "Ddit3","Trib3") #2nd Adjustment Sophie October 2019 (V3)

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Amigo_chol_genes_known, showlabels = F, Coi=c("")) + 
  ggtitle("Amigo cholesterol genes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Amigo_chol_genes_known, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Amigo_chol_genes_known, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_Amigo_chol_genes_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Amigo_chol_genes_known, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_Amigo_chol_genes_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(Amigo_chol_genes_known, Gdiffexp)
nonSigGenes<-setdiff(Amigo_chol_genes_known, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Amigo_chol_genes_DE_3", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Amigo_chol_genes_nonDE_3", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Amigo_chol_genes_nonExp_3", append = T)
}

########################################
##### DC function
########################################

# ##### DC maturation ====> DECIDED TO REMOVE
# DCmaturationGenes<-c("ATF2","CCR7","Cd1","CD1A","CD40","CD40LG","CD58","CD80","CD83","CD86","Creb","CSF2","DDR2","ERK1","ERK2","Fascin",
#                      "Fcgr1","Fcgr2","FCGR2B","HLA-DR","ICAM1","IFNa","IFNb","IFNAR1","IgG","Ikb","IkB-NfkB","IKK","IL1","IL10","IL12",
#                      "IL12B","IL15","IL1B","IL1RL2","IL23A","IL32","IL36B","IL36G","IL6")
# DCmaturationGenes<-firstup(tolower(DCmaturationGenes))
# setdiff(DCmaturationGenes, rownames(expTable_mean))
# setdiff(DCmaturationGenes, rownames(countData))
# DCmaturationGenes[DCmaturationGenes=="Cd1a"]<-"Cd1d1"
# DCmaturationGenes[DCmaturationGenes=="Creb"]<-"Creb1"
# DCmaturationGenes[DCmaturationGenes=="Erk1"]<-"Mapk3"
# DCmaturationGenes[DCmaturationGenes=="Erk2"]<-"Mapk1"
# DCmaturationGenes[DCmaturationGenes=="Fascin"]<-"Fscn1"
# DCmaturationGenes[DCmaturationGenes=="Ifna"]<-"Ifna1"
# DCmaturationGenes<-c(DCmaturationGenes,"Ifna2")
# DCmaturationGenes<-c(DCmaturationGenes,"Ifna15")
# DCmaturationGenes[DCmaturationGenes=="Ifnb"]<-"Ifnb1"
# DCmaturationGenes[DCmaturationGenes=="Ikb"]<-"Nfkbib"
# DCmaturationGenes[DCmaturationGenes=="Ikb-nfkb"]<-"Nfkbia"
# DCmaturationGenes[DCmaturationGenes=="Ikk"]<-"Ikbkb"
# DCmaturationGenes<-c(DCmaturationGenes,"Ikbkg")
# DCmaturationGenes[DCmaturationGenes=="Il1"]<-"Il1a"
# DCmaturationGenes<-c(DCmaturationGenes,"Il1b")
# DCmaturationGenes<-c(DCmaturationGenes,"Il1rn")
# DCmaturationGenes[DCmaturationGenes=="Il12"]<-"Il12a"
# DCmaturationGenes<-c(DCmaturationGenes,"Il12b")
# DCmaturationGenes[DCmaturationGenes=="Il36b"]<-"Il1f8"
# DCmaturationGenes[DCmaturationGenes=="Il36g"]<-"Il1f9"
# DCmaturationGenes<-intersect(DCmaturationGenes, rownames(expTable_mean))
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationGenes, showlabels = F, Coi=c("")) + 
#   ggtitle("DCmaturationGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationGenes, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)

##### DC maturation, TLR induced
genesC49<-c("CMPK2","OAS3","ISG20","SLFN8","NT5C3","IFI47","IFIT2","DDX58","SLFN3","IFIT1","SLFN1","I830012O16RIK","IRF7","OASL2","SCO1",
            "TOR3A","PML","SLFN5","DAXX","EIF2AK2","DHX58","H2-T24","OAS1A","RTP4","PARP12","TMEM184B","SLC25A22","1600014C10RIK","XAF1",
            "B4GALT5","IFI44","TRIM30A","IFI204","MTHFR","STAT2","RNF213","LGALS9","ATRIP","IRGM1","LGALS3BP","5730508B09RIK","SGCB",
            "TOR1AIP1","MS4A4C","CSPRS","DTX3L","ZBP1","PDE7B","AI607873","USP18","ACADL","GVIN1","SP100","LTBP1","CARHSP1","GMPPB","KDR",
            "IRF9","IFI203","LY6A","LY6C2","RNF4","PPM1K","DR1","SLFN2","CDH1","ZBTB5","CDC42BPG","CHIC1","AI837181","VCPIP1","CCDC23",
            "TRIM25","ORAI1","SLC43A3","MELK","4632434I11RIK")
genesC61<-c("ARID5A","BBX","RHBDF2","IFIH1","PARP14","OASL1","RAPGEF2","SPHK1","ABTB2","CD69","CASP4","PARP9","RSAD2","XRN1","FABP7","RIPK2",
            "AVL9","OAS1B","ZUFSP","CMTR1","TDRD7","GBP7","PARP10","ZNFX1","CCDC6","PIK3R5","BATF2","TOR1AIP2","MIR155","NAMPT","SLFN10-PS",
            "CH25H","FBXO42","IL6","OAS2","IL2RA","SMG7","NR4A3","MUL1","ADAR","FNDC3A","CDKN1A","PCSK7","CXCL11","PLAGL1","ADAP2","DNAJC13",
            "D730005E14RIK","SETDB2","FAP","FAM84B","ASTL","NDRG1","MX1","LRRC63","BC016423","TJAP1","CD80","AA467197","WHSC1L1","IRG1",
            "NUP210L","TLCD1","GPR85","EIF1A","MUC1","IL10","KIF1A","SMCR8","ZFP318","CREM","SLC30A4","MID1","5031425E22RIK","ATP6V1A")
DCmaturationTLR<-firstup(tolower(c(genesC49, genesC61)))
setdiff(DCmaturationTLR, rownames(expTable_mean))
setdiff(DCmaturationTLR, rownames(countData))
DCmaturationTLR[DCmaturationTLR=="I830012o16rik"]<-"Ifit3b"
DCmaturationTLR[DCmaturationTLR=="H2-t24"]<-"H2-T24"
DCmaturationTLR[DCmaturationTLR=="1600014c10rik"]<-"1600014C10Rik"
DCmaturationTLR[DCmaturationTLR=="5730508b09rik"]<-"Fam241a"
DCmaturationTLR[DCmaturationTLR=="Ai607873"]<-"Ifi207"
DCmaturationTLR[DCmaturationTLR=="Ai837181"]<-"AI837181"
DCmaturationTLR[DCmaturationTLR=="4632434i11rik"]<-"Ddias"
DCmaturationTLR[DCmaturationTLR=="D730005e14rik"]<-"D730005E14Rik"
DCmaturationTLR[DCmaturationTLR=="Bc016423"]<-"Fam208b"
DCmaturationTLR[DCmaturationTLR=="Aa467197"]<-"AA467197"
DCmaturationTLR[DCmaturationTLR=="5031425e22rik"]<-"5031425E22Rik"
nonExpGenes <-setdiff(DCmaturationTLR, rownames(expTable_mean))
nonExpGenes2 <-setdiff(DCmaturationTLR, rownames(expTable_mean))
DCmaturationTLR<-intersect(DCmaturationTLR, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationTLR, showlabels = F, Coi=c("")) + 
  ggtitle("DCmaturationTLR") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationTLR, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationTLR, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationTLR_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationTLR, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationTLR_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationTLR, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationTLR, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationTLR_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationTLR_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationTLR_nonExp", append = T)
}

##### DC maturation, homeostatic
genesC103<-c("SLCO5A1","RASAL2","GALNT12","TBC1D15","MMP7","MICAL3","1700021K19RIK","SDHAF1","MMP23","HLX","NUAK1","ENO3","PIK3R3",
             "JOSD1","VHL","HECW2","KCNN1","CHKA","SCUBE3","ZFP36L1","IFT57","ETFDH","APPL1","DNAJC22","CSNK1G1","CABLES2","CHST2",
             "2200002D01RIK","TPM1","LZTS2","MPP5","ACE2","SQLE","CDCP1","TSC1","FUT4","CRELD1","KIF13B","ACP2","WDR91","NOSTRIN",
             "SPATS2","B4GALNT1","TMEM8","LEPRE1","TDH","GFPT1","ST8SIA6","CCDC120","ELK3","ARFRP1","STXBP2","HMGCS1","A430093F15RIK",
             "RAB30","RAD50","TNNT2","CYP51","SC5D","PHKA1","CLN5","ZFP180","AK2","MYLIP","VCAM1","SLC7A6","RCOR3","ZFP568","TTC39A",
             "MVD","FBXL3","GLT28D2","SOX4","PYROXD1","NUP62","KDM5B","INCA1","PRDX3","MKKS","VRK2","MIOS","PPCS","CD2BP2","RB1CC1",
             "1810033B17RIK","PSMD1","AI314180","PSCA")
DCmaturationHomeo<-firstup(tolower(genesC103))
DCmaturationHomeo<-c(DCmaturationHomeo,"Nuak1","Pik3r3","Hecw2","Kcnn1","Scube3","Cdcp1","Ccdc120") #Adjustment Sophie
setdiff(DCmaturationHomeo, rownames(expTable_mean))
setdiff(DCmaturationHomeo, rownames(countData))
DCmaturationHomeo[DCmaturationHomeo=="1700021k19rik"]<-"Rubcn"
DCmaturationHomeo[DCmaturationHomeo=="2200002d01rik"]<-"2200002D01Rik"
DCmaturationHomeo[DCmaturationHomeo=="Lepre1"]<-"P3h1"
DCmaturationHomeo[DCmaturationHomeo=="A430093f15rik"]<-"A430093F15Rik"
DCmaturationHomeo[DCmaturationHomeo=="1810033b17rik"]<-"Mcemp1"
DCmaturationHomeo[DCmaturationHomeo=="Ai314180"]<-"AI314180"
nonExpGenes <-setdiff(DCmaturationHomeo, rownames(expTable_mean))
nonExpGenes2 <-c(nonExpGenes2, setdiff(DCmaturationHomeo, rownames(expTable_mean))) #For merged list
DCmaturationHomeo<-intersect(DCmaturationHomeo, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationHomeo, showlabels = F, Coi=c("")) + 
  ggtitle("DCmaturationHomeo") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationHomeo, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationHomeo, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationHomeo_test_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationHomeo, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationHomeo_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationHomeo, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationHomeo, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationHomeo_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationHomeo_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationHomeo_nonExp", append = T)
}

##### DC maturation, common immunogenic + homeostatic
genesC72<-c("RNF31","NIPAL1","GRAMD3","RAB22A","KPNA3","PIK3R1","MFHAS1","PDCD1LG2","MOCOS","GYPC","2810474O19RIK","DRAM1","GBP2","GLIPR2",
            "RNF24","ZFP513","KDM5C","MAPKAPK3","LPP","CD274","STAT3","HERC6","TMEM106A","MIER3","ANKIB1","EHD1","ZFP281","CCL5","NFKBIB",
            "ARHGAP8","CAR13","RELA","CD40","GBP3","AKT3","1110032F04RIK","TBL1X","CSRP1","CASZ1","TAP1","LRP8","PLA1A","INPP5B","WSB2",
            "ELMOD2","FZD5","NUPR1","AGRN","SUSD2","OSM","HSH2D","LITAF","ANXA4","RASGRP1","TLE3","RECK","ETNK1","POLD3","OPTN","BCL3",
            "PCGF3","UBE2L6","STAT1","STAT4","GRINA","2610002M06RIK","NINJ1","SOCS1","PIAS1","SERPINB1B","OCLN","PLAGL2","CXCR5","BCL2L11",
            "AIDA","BC002230","C130026I21RIK","SNX10","SFMBT2","RAB8B","CAR2","CNOT4","RND1","TRAFD1","NDRG2","FLNB","SERTAD3",
            "9330175E14RIK","HSD17B11","PYHIN1","RUSC1","LAMP2","ISPD","DAPP1","GGTA1","VDR","TAS1R3","PHF16","CNN3","BC016579","CDC42EP4",
            "S1PR3","TAB2","ZBTB22")
genesC73<-c("NFKB2","SPRED1","SIK3","AEBP2","NUAK2","DIRC2","DUSP5","SLC2A6","TNIP1","MIF4GD","LNX2","STX11","KHNYN","SAV1","PVR","SERINC5",
            "MXD1","FAM126A","ATP11A","PRNP","GPR126","GADD45B","MVP","ETV3","BIRC3","1700047I17RIK2","TRIP10","PTGER4","UBXN2A","TMCC3",
            "BASP1","GABPB1","PHLPP1","CEP350","ARHGAP31","PHLDB1","TRAF6","ATP2B1","RAB12","RABGEF1","MYO1C","RNF19B","AFTPH","ATP6V0A1",
            "PI4K2B","CERS6","POGK","SIN3B","FNBP1L","PPP4R2","MLLT6","LY75","SEMA4C","MTMR12","CTTNBP2NL","NFKBIA","ICAM1","SLC6A6",
            "ARHGEF3","GM527","FRMD4A","LRRK1","SEMA7A","STXBP3A","DENND3","MAP2K1","STAT5A","MARCKSL1","TMTC2","DCBLD2","SBNO2","BCL2L1",
            "PLEKHG2","AP1G2","CLIC4","MT2","SLC22A15","LRIG2","LLGL1","SMARCE1","LDLR","INPP5A","STK19","HERC3","OTUD5","SMG1","CCL22",
            "NBL1","WDR37","D17WSU92E","TMCO4","KDM2B","NRP2","SLC33A1","SDC4","RASA4","TAB3","EXOC6B","HBP1","TNFAIP2","PRDM1","MAPKBP1",
            "H2-T10","TANC1","IL12B","MFSD1","MAGI3","NCOA5","NPEPPS","LGALS8","PPP3CB")
genesC115<-c("IL15RA","TMEM39A","C030046E11RIK","ARHGAP28","SWAP70","BIRC2","FAS","SNN","INSL6","TNFRSF4","CLEC2I","PCGF5","IL2RG","LOXL3",
             "IFFO2","MDM2","ZDHHC14","PAPSS2","ANKRD33B","IL4RA","C230081A13RIK","GTDC1","PEX13","IL21R","ADORA2A","TANK","NCOA7","AMN1",
             "CD200","UAP1","TRAF1","GBP8","PAQR3","ATXN1","ZBTB10","APOBEC1","IL4I1","TAPBPL","THAP2","SPG11","TCF7","RHEBL1","CCR7","SPECC1",
             "NSUN4","SNX11","TMEM176A","ARHGAP22","CBLB","SIRT1","LACTB","CDC42EP3","IL15","PRR14","BCL2L14","ASPRV1","NUDT9","RRAD","ZBTB37",
             "FAM53B","ST6GALNAC6","TMTC4","NABP1","TBC1D1","D930015E06RIK","RNF43","PURG","3110043O21RIK","NXF1","TMEM140","MTF1","KCNK6",
             "LAMP1","DOCK9","FOXP4","KTN1","RRAS2","GCA","NFKBIE","SPIB","SNX16","PDZK1IP1","N4BP2L1","LMBR1L","ZHX2","PARP3","HOOK2",
             "B3GNT5","FABP5","INO80D","ATXN7L1","SERPINB9","OSBPL8","PAPOLG","EHF","KRIT1","GPR132","SHMT1","HSF2","CFLAR","EPS8L2","WRN",
             "MFSD7B","RENBP","ANKRD12","ART2B","CPNE2","IRF1","RASIP1","INPP5F","DUSP14","TMEM70","PHXR4","DAO","KATNA1","CLEC2G","MAP1B",
             "MED11","CASP1","5730419I09RIK","RAB20","PDE2A","CLIP1")
DCmaturationCommon<-firstup(tolower(c(genesC72, genesC73, genesC115)))
setdiff(DCmaturationCommon, rownames(expTable_mean))
setdiff(DCmaturationCommon, rownames(countData))
DCmaturationCommon[DCmaturationCommon=="2810474o19rik"]<-"2810474O19Rik"
DCmaturationCommon[DCmaturationCommon=="1110032f04rik"]<-"1110032F04Rik"
DCmaturationCommon[DCmaturationCommon=="2610002m06rik"]<-"2610002M06Rik"
DCmaturationCommon[DCmaturationCommon=="Bc002230"]<-"Nrde2"
DCmaturationCommon[DCmaturationCommon=="C130026i21rik"]<-"C130026I21Rik"
DCmaturationCommon[DCmaturationCommon=="9330175e14rik"]<-"9330175E14Rik"
DCmaturationCommon[DCmaturationCommon=="Phf16"]<-"Jade3"
DCmaturationCommon[DCmaturationCommon=="Bc016579"]<-"BC016579"
DCmaturationCommon[DCmaturationCommon=="Gpr126"]<-"Adgrg6"
DCmaturationCommon[DCmaturationCommon=="1700047i17rik2"]<-"1700047I17Rik2"
DCmaturationCommon[DCmaturationCommon=="Stxbp3a"]<-"Nsfl1c"
DCmaturationCommon[DCmaturationCommon=="D17wsu92e"]<-"D17Wsu92e"
DCmaturationCommon[DCmaturationCommon=="H2-t10"]<-"H2-T10"
DCmaturationCommon[DCmaturationCommon=="C030046e11rik"]<-"Ric1"
DCmaturationCommon[DCmaturationCommon=="C230081a13rik"]<-"Peak1"
DCmaturationCommon[DCmaturationCommon=="D930015e06rik"]<-"Tmem131l"
DCmaturationCommon[DCmaturationCommon=="3110043o21rik"]<-"3110043O21Rik"
DCmaturationCommon[DCmaturationCommon=="5730419i09rik"]<-"C2cd5"
nonExpGenes <-setdiff(DCmaturationCommon, rownames(expTable_mean))
nonExpGenes2 <-c(nonExpGenes2, setdiff(DCmaturationCommon, rownames(expTable_mean))) #For merged list
DCmaturationCommon<-intersect(DCmaturationCommon, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommon, showlabels = F, Coi=c("")) + 
  ggtitle("DCmaturationCommon") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommon, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommon, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationCommon_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommon, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationCommon_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationCommon, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationCommon, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationCommon_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationCommon_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationCommon_nonExp", append = T)
}

##### DC maturation, common immunogenic + homeostatic (new: January 2020)
genesC72<-c("RNF31","NIPAL1","GRAMD3","RAB22A","KPNA3","PIK3R1","MFHAS1","PDCD1LG2","MOCOS","GYPC","2810474O19RIK","DRAM1","GBP2","GLIPR2",
            "RNF24","ZFP513","KDM5C","MAPKAPK3","LPP","CD274","STAT3","HERC6","TMEM106A","MIER3","ANKIB1","EHD1","ZFP281","CCL5","NFKBIB",
            "ARHGAP8","CAR13","RELA","CD40","GBP3","AKT3","1110032F04RIK","TBL1X","CSRP1","CASZ1","TAP1","LRP8","PLA1A","INPP5B","WSB2",
            "ELMOD2","FZD5","NUPR1","AGRN","SUSD2","OSM","HSH2D","LITAF","ANXA4","RASGRP1","TLE3","RECK","ETNK1","POLD3","OPTN","BCL3",
            "PCGF3","UBE2L6","STAT1","STAT4","GRINA","2610002M06RIK","NINJ1","SOCS1","PIAS1","SERPINB1B","OCLN","PLAGL2","CXCR5","BCL2L11",
            "AIDA","BC002230","C130026I21RIK","SNX10","SFMBT2","RAB8B","CAR2","CNOT4","RND1","TRAFD1","NDRG2","FLNB","SERTAD3",
            "9330175E14RIK","HSD17B11","PYHIN1","RUSC1","LAMP2","ISPD","DAPP1","GGTA1","VDR","TAS1R3","PHF16","CNN3","BC016579","CDC42EP4",
            "S1PR3","TAB2","ZBTB22")
genesC73<-c("NFKB2","SPRED1","SIK3","AEBP2","NUAK2","DIRC2","DUSP5","SLC2A6","TNIP1","MIF4GD","LNX2","STX11","KHNYN","SAV1","PVR","SERINC5",
            "MXD1","FAM126A","ATP11A","PRNP","GPR126","GADD45B","MVP","ETV3","BIRC3","1700047I17RIK2","TRIP10","PTGER4","UBXN2A","TMCC3",
            "BASP1","GABPB1","PHLPP1","CEP350","ARHGAP31","PHLDB1","TRAF6","ATP2B1","RAB12","RABGEF1","MYO1C","RNF19B","AFTPH","ATP6V0A1",
            "PI4K2B","CERS6","POGK","SIN3B","FNBP1L","PPP4R2","MLLT6","LY75","SEMA4C","MTMR12","CTTNBP2NL","NFKBIA","ICAM1","SLC6A6",
            "ARHGEF3","GM527","FRMD4A","LRRK1","SEMA7A","STXBP3A","DENND3","MAP2K1","STAT5A","MARCKSL1","TMTC2","DCBLD2","SBNO2","BCL2L1",
            "PLEKHG2","AP1G2","CLIC4","MT2","SLC22A15","LRIG2","LLGL1","SMARCE1","LDLR","INPP5A","STK19","HERC3","OTUD5","SMG1","CCL22",
            "NBL1","WDR37","D17WSU92E","TMCO4","KDM2B","NRP2","SLC33A1","SDC4","RASA4","TAB3","EXOC6B","HBP1","TNFAIP2","PRDM1","MAPKBP1",
            "H2-T10","TANC1","IL12B","MFSD1","MAGI3","NCOA5","NPEPPS","LGALS8","PPP3CB")
genesC115<-c("IL15RA","TMEM39A","C030046E11RIK","ARHGAP28","SWAP70","BIRC2","FAS","SNN","INSL6","TNFRSF4","CLEC2I","PCGF5","IL2RG","LOXL3",
             "IFFO2","MDM2","ZDHHC14","PAPSS2","ANKRD33B","IL4RA","C230081A13RIK","GTDC1","PEX13","IL21R","ADORA2A","TANK","NCOA7","AMN1",
             "CD200","UAP1","TRAF1","GBP8","PAQR3","ATXN1","ZBTB10","APOBEC1","IL4I1","TAPBPL","THAP2","SPG11","TCF7","RHEBL1","CCR7","SPECC1",
             "NSUN4","SNX11","TMEM176A","ARHGAP22","CBLB","SIRT1","LACTB","CDC42EP3","IL15","PRR14","BCL2L14","ASPRV1","NUDT9","RRAD","ZBTB37",
             "FAM53B","ST6GALNAC6","TMTC4","NABP1","TBC1D1","D930015E06RIK","RNF43","PURG","3110043O21RIK","NXF1","TMEM140","MTF1","KCNK6",
             "LAMP1","DOCK9","FOXP4","KTN1","RRAS2","GCA","NFKBIE","SPIB","SNX16","PDZK1IP1","N4BP2L1","LMBR1L","ZHX2","PARP3","HOOK2",
             "B3GNT5","FABP5","INO80D","ATXN7L1","SERPINB9","OSBPL8","PAPOLG","EHF","KRIT1","GPR132","SHMT1","HSF2","CFLAR","EPS8L2","WRN",
             "MFSD7B","RENBP","ANKRD12","ART2B","CPNE2","IRF1","RASIP1","INPP5F","DUSP14","TMEM70","PHXR4","DAO","KATNA1","CLEC2G","MAP1B",
             "MED11","CASP1","5730419I09RIK","RAB20","PDE2A","CLIP1")
genesc114 <- c("VSIG10","CCSER2","TMEM150C","ENO2","RCSD1","GATSL2","TMEM19","JAG1","RELB","EXTL1",
               "STXBP1","TSPAN3","SSH1","SH3BP4","DAAM1","TNFRSF11A","AGAP1", "ZC3H7B","ARC","GAL3ST2",
               "E130311K13RIK","POGLUT1","ADCY6", "PIAS3", "STK4", "RNF115","SAMSN1","SLC26A2","RAP2B", "NUDT17",
               "PLCL2", "ATRNL1","LIMA1", "IDI1","BMP2K","CACNB3","MEX3B", "TRAF2", "PFKFB3","KDM1B",
               "ANXA3", "STOML1","ZFAND6","CDK17", "RASA2","CDKN2B","RAB21", "CREBL2","ENOX2", "TNFRSF1B",
               "CCNG2", "NFE2L1","SOCS2", "IL7R","DOK1","SLC22A23","SLC4A8","MYO1G", "TMEM123","HCN3", 
               "EFNA2", "PPP1R16A","PNPLA8","SPSB1", "LAD1","RNF19A","GRAMD1B","POLR3C","FAH","STARD7",
               "SPINT2","ZMAT1", "RFTN1", "PTPN21","RBPMS","TNFRSF9","STAP2", "TMEM176B","MREG","ADAP1",
               "SEL1L", "4632411B12RIK","SERPINB9B","PLEKHM2","INF2","GCLC","CLPTM1","ATMIN", "ASL","RASSF3",
               "ZBTB18","ZFP46", "C330006D17RIK","SCIN","NFAT5","RABGAP1L","LLGL2", "DIP2B", "HERC4", "PLA2G4F",
               "MAP3K14","RGS12", "VPS13D","AKIRIN2","FSCN1","KSR1","KIF1B", "IDH1","GNB4","RALGAPA2",
               "CSF1","SIDT2", "ARHGEF40","CDK13", "DNM1L","TMEM18","SEPSECS","RSPRY1","AP3M2", "ADM",
               "ZFYVE1","SEMA6D","DONSON","1700081L11RIK","CSF2RB","H2-EB2","ASPSCR1","FAM129C","TBC1D10A","NCF1", 
               "TADA1", "EML6","KCTD9", "FBXO7", "CHST7","DCLRE1A","ULK3","ANKRD27","CSF2RB2","ZW10", 
               "CDC14B","DDHD1", "NUP85", "FOXH1", "IKBKE","NR2F6", "SLAMF1","TRAF5", "MTMR4", "FBRSL1",
               "WBSCR27","NCK1","STYX","S1PR1", "GET4","NEDD1")

DCmaturationCommonNew<-firstup(tolower(c(genesC72, genesC73, genesC115,genesc114)))
setdiff(DCmaturationCommonNew, rownames(expTable_mean))
setdiff(DCmaturationCommonNew, rownames(countData))
DCmaturationCommonNew[DCmaturationCommonNew=="2810474o19rik"]<-"2810474O19Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="1110032f04rik"]<-"1110032F04Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="2610002m06rik"]<-"2610002M06Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="Bc002230"]<-"Nrde2"
DCmaturationCommonNew[DCmaturationCommonNew=="C130026i21rik"]<-"C130026I21Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="9330175e14rik"]<-"9330175E14Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="Phf16"]<-"Jade3"
DCmaturationCommonNew[DCmaturationCommonNew=="Bc016579"]<-"BC016579"
DCmaturationCommonNew[DCmaturationCommonNew=="Gpr126"]<-"Adgrg6"
DCmaturationCommonNew[DCmaturationCommonNew=="1700047i17rik2"]<-"1700047I17Rik2"
DCmaturationCommonNew[DCmaturationCommonNew=="Stxbp3a"]<-"Nsfl1c"
DCmaturationCommonNew[DCmaturationCommonNew=="D17wsu92e"]<-"D17Wsu92e"
DCmaturationCommonNew[DCmaturationCommonNew=="H2-t10"]<-"H2-T10"
DCmaturationCommonNew[DCmaturationCommonNew=="C030046e11rik"]<-"Ric1"
DCmaturationCommonNew[DCmaturationCommonNew=="C230081a13rik"]<-"Peak1"
DCmaturationCommonNew[DCmaturationCommonNew=="D930015e06rik"]<-"Tmem131l"
DCmaturationCommonNew[DCmaturationCommonNew=="3110043o21rik"]<-"3110043O21Rik"
DCmaturationCommonNew[DCmaturationCommonNew=="5730419i09rik"]<-"C2cd5"
DCmaturationCommonNew[DCmaturationCommonNew=="4632411b12rik"]<-"Kansl3"
DCmaturationCommonNew[DCmaturationCommonNew=="C330006d17rik"]<-"Birc2"
DCmaturationCommonNew[DCmaturationCommonNew=="1700081l11rik"]<-"Kansl1"
nonExpGenes <-setdiff(DCmaturationCommonNew, rownames(expTable_mean))
nonExpGenes2 <-c(nonExpGenes2, setdiff(DCmaturationCommonNew, rownames(expTable_mean))) #For merged list
DCmaturationCommonNew<-intersect(DCmaturationCommonNew, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F, Coi=c("")) + 
  ggtitle("DCmaturationCommon") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationCommonNew_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationCommonNew_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationCommonNew, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationCommonNew, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/Suppl_Table3_summary_ownGOlists_final_paper_jan2020_formatted_V2.xlsx", "DCmaturationCommonNew_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/Suppl_Table3_summary_ownGOlists_final_paper_jan2020_formatted_V2.xlsx", "DCmaturationCommonNew_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/Suppl_Table3_summary_ownGOlists_final_paper_jan2020_formatted_V2.xlsx", "DCmaturationCommonNew_nonExp", append = T)
}

##### DC maturation extra
genesc114 <- read.xlsx("../documentation/DC_maturation_data.xlsx", sheetName = "C114")
genesc114 <- genesc114$Gene.Symbol
as.character(genesc114)
genesc114 <- c("VSIG10","CCSER2","TMEM150C","ENO2","RCSD1","GATSL2","TMEM19","JAG1","RELB","EXTL1",
"STXBP1","TSPAN3","SSH1","SH3BP4","DAAM1","TNFRSF11A","AGAP1", "ZC3H7B","ARC","GAL3ST2",
"E130311K13RIK","POGLUT1","ADCY6", "PIAS3", "STK4", "RNF115","SAMSN1","SLC26A2","RAP2B", "NUDT17",
"PLCL2", "ATRNL1","LIMA1", "IDI1","BMP2K","CACNB3","MEX3B", "TRAF2", "PFKFB3","KDM1B",
"ANXA3", "STOML1","ZFAND6","CDK17", "RASA2","CDKN2B","RAB21", "CREBL2","ENOX2", "TNFRSF1B",
"CCNG2", "NFE2L1","SOCS2", "IL7R","DOK1","SLC22A23","SLC4A8","MYO1G", "TMEM123","HCN3", 
"EFNA2", "PPP1R16A","PNPLA8","SPSB1", "LAD1","RNF19A","GRAMD1B","POLR3C","FAH","STARD7",
"SPINT2","ZMAT1", "RFTN1", "PTPN21","RBPMS","TNFRSF9","STAP2", "TMEM176B","MREG","ADAP1",
"SEL1L", "4632411B12RIK","SERPINB9B","PLEKHM2","INF2","GCLC","CLPTM1","ATMIN", "ASL","RASSF3",
"ZBTB18","ZFP46", "C330006D17RIK","SCIN","NFAT5","RABGAP1L","LLGL2", "DIP2B", "HERC4", "PLA2G4F",
"MAP3K14","RGS12", "VPS13D","AKIRIN2","FSCN1","KSR1","KIF1B", "IDH1","GNB4","RALGAPA2",
"CSF1","SIDT2", "ARHGEF40","CDK13", "DNM1L","TMEM18","SEPSECS","RSPRY1","AP3M2", "ADM",
"ZFYVE1","SEMA6D","DONSON","1700081L11RIK","CSF2RB","H2-EB2","ASPSCR1","FAM129C","TBC1D10A","NCF1", 
"TADA1", "EML6","KCTD9", "FBXO7", "CHST7","DCLRE1A","ULK3","ANKRD27","CSF2RB2","ZW10", 
"CDC14B","DDHD1", "NUP85", "FOXH1", "IKBKE","NR2F6", "SLAMF1","TRAF5", "MTMR4", "FBRSL1",
"WBSCR27","NCK1","STYX","S1PR1", "GET4","NEDD1")
DCmaturationExtra<-firstup(tolower(genesc114))
setdiff(DCmaturationExtra, rownames(expTable_mean))
setdiff(DCmaturationExtra, rownames(countData))
DCmaturationExtra[DCmaturationExtra=="4632411b12rik"]<-"Kansl3"
DCmaturationExtra[DCmaturationExtra=="C330006d17rik"]<-"Birc2"
DCmaturationExtra[DCmaturationExtra=="1700081l11rik"]<-"Kansl1"

nonExpGenes <-setdiff(DCmaturationExtra, rownames(expTable_mean))
nonExpGenes2 <-c(nonExpGenes2, setdiff(DCmaturationExtra, rownames(expTable_mean))) #For merged list
DCmaturationExtra<-intersect(DCmaturationExtra, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationExtra, showlabels = F, Coi=c("")) + 
  ggtitle("DCmaturationExtra") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationExtra, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationExtra, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationExtra_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationExtra, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationExtra_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationExtra, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationExtra, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationExtra_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationExtra_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationExtra_nonExp", append = T)
}


##### DC maturation merge 
DCmaturationMerge<-unique(c(DCmaturationCommon,DCmaturationExtra,DCmaturationHomeo,DCmaturationTLR))
setdiff(DCmaturationMerge, rownames(expTable_mean))
setdiff(DCmaturationMerge, rownames(countData))
nonExpGenes2 #Check merged list
DCmaturationMerge<-intersect(DCmaturationMerge, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationMerge, showlabels = F, Coi=c("")) +
  ggtitle("DCmaturationMerge") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationMerge, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationMerge, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_DCmaturationMerge_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationMerge, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_DCmaturationMerge_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(DCmaturationMerge, Gdiffexp)
nonSigGenes<-setdiff(DCmaturationMerge, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationMerge_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationMerge_nonDE", append = T)
if (length(nonExpGenes2) > 1) {
  write.xlsx(nonExpGenes2, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "DCmaturationMerge", append = T)
}



# ##### Negative regulation inflammation
# inflammationGenes<-unique(c(gsets$`GO:001818`,gsets$`GO:50728`))
# setdiff(inflammationGenes, rownames(expTable_mean))
# setdiff(inflammationGenes, rownames(countData))
# inflammationGenes<-intersect(inflammationGenes, rownames(expTable_mean))
# 
# ## GO:001818 = ??
# ## GO:50728 = ??
# 
# p1<-plotDotplot(barycoords, Gdiffexp, Goi=inflammationGenes, showlabels = F, Coi=c("")) + 
#   ggtitle("inflammationGenes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# p2<-plotRoseplot(barycoords, Gdiffexp, Goi=inflammationGenes, showlabels = F, Coi=c(""))
# grid.arrange(p1, p2, ncol=2)


########################################
##### Apoptotic engulfment receptors
########################################
# engulfmentGenes<-c("Ucp2","Lrrn1","Mfge8","Adgrf3","Adgrf5","Adora2a","Clec2i","Slco5a1","Src","C1qtnf6","Slco2a1","Megf11",
#                    "Adgrb1","Timd4","Stab1","Stab2","Ager","Cd300lf","Trem2","Mertk","Tyro3","Axl","Gas6","Cd36","Lrp1","C1qA",
#                    "Scarf","Bai3","Bai1","Havcr","Avb5","Avb3","Abca1")
engulfmentGenes<-c("Slco5a1","Src","Clec2i","Adgrf5","Mfge8","Adora2a","Ucp2")
setdiff(engulfmentGenes, rownames(expTable_mean))
setdiff(engulfmentGenes, rownames(countData))
# engulfmentGenes[engulfmentGenes=="C1qA"]<-"C1qa"
# engulfmentGenes[engulfmentGenes=="Bai3"]<-"Adgrb3"
# engulfmentGenes[engulfmentGenes=="Bai1"]<-"Adgrb1"
# engulfmentGenes[engulfmentGenes=="Havcr"]<-"Havcr1"

engulfmentGenes<-unique(engulfmentGenes)
nonExpGenes <-setdiff(engulfmentGenes, rownames(expTable_mean))
engulfmentGenes<-intersect(engulfmentGenes, rownames(expTable_mean))

##### Added on 11/09/2018 #####
engulfmentGenes<-c("Ucp2","Lrrn1","Mfge8","Adgrf3","Adgrf5","Adora2a","Clec2i","Slco5a1","Src","C1qtnf6","Slco2a1","Megf11","Adgrb1",
                   "Timd4","Stab1","Stab2","Ager","Cd300lf","Trem2","Mertk","Tyro3","Axl","Gas6","Cd36","Lrp1","C1qA","Scarf","Bai3",
                   "Bai1","Havcr","Avb5","Avb3","Abca1")
engulfmentGenes<-c(engulfmentGenes,"Itgae","Cd207") ##extra addition 23/09/19 at request Victor and Sophie
setdiff(engulfmentGenes, rownames(expTable_mean))
setdiff(engulfmentGenes, rownames(countData))
engulfmentGenes[engulfmentGenes=="C1qA"]<-"C1qa"
engulfmentGenes[engulfmentGenes=="Scarf"]<-"Calm4"
engulfmentGenes[engulfmentGenes=="Bai3"]<-"Adgrb3"
engulfmentGenes[engulfmentGenes=="Bai1"]<-"Adgrb1"
engulfmentGenes[engulfmentGenes=="Havcr"]<-"Havcr1"

engulfmentGenes<-unique(engulfmentGenes)
nonExpGenes <-setdiff(engulfmentGenes, rownames(expTable_mean))
engulfmentGenes<-intersect(engulfmentGenes, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=engulfmentGenes, showlabels = F, Coi=c("")) + 
  ggtitle("Apoptopic engulment") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=engulfmentGenes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=engulfmentGenes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_apoptopicEngulfment_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=engulfmentGenes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_apoptopicEngulfment_final_paper.png",dpi=300)

### Split in DE and non-DE genes
sigGenes<-intersect(engulfmentGenes, Gdiffexp)
nonSigGenes<-setdiff(engulfmentGenes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "engulfmentGenes_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "engulfmentGenes_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "engulfmentGenes_nonExp", append = T)
}



########################################
##### Apoptotic Cell Clearance (Clint)
########################################

#Read in Amigo results from excel file
Apop_cell_clearance <- read_xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/documentation/Direct_GO_Apoptotic_cell_clearance.xlsx", col_names = F)
Apop_cell_clearance_full <- read_xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/documentation/Full_GO_Apoptotic_cell_clearance.xlsx", col_names = F)

Apop_cell_clearance_genes<-unique(Apop_cell_clearance[[1]])
Apop_cell_clearance_genes_full<-unique(Apop_cell_clearance_full[[1]])


setdiff(Apop_cell_clearance_genes, rownames(expTable_mean))
setdiff(Apop_cell_clearance_genes, rownames(countData))

nonExpGenes <-setdiff(Apop_cell_clearance_genes, rownames(expTable_mean))
Apop_cell_clearance_genes<-intersect(Apop_cell_clearance_genes, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Apop_cell_clearance_genes, showlabels = F, Coi=c("")) + 
  ggtitle("Apoptotic Clearance genes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Apop_cell_clearance_genes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Apop_cell_clearance_genes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_Apoptotic_Clearance_genes_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Apop_cell_clearance_genes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_Apoptotic_Clearance_genes_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(Apop_cell_clearance_genes, Gdiffexp)
nonSigGenes<-setdiff(Apop_cell_clearance_genes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_nonExp", append = T)
}


###########################################
##### Engulfment of Apoptotic Cell  (Clint)
###########################################

#Read in Amigo results from excel file
Engulfment_clint <- read_xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/documentation/Direct_GO_Engulfment_of_Apoptotic_cell.xlsx", col_names = F)

Engulfment_genes_clint<-unique(Engulfment_clint[[1]])

setdiff(Engulfment_genes_clint, rownames(expTable_mean))
setdiff(Engulfment_genes_clint, rownames(countData))

nonExpGenes <-setdiff(Engulfment_genes_clint, rownames(expTable_mean))
Engulfment_genes_clint<-intersect(Engulfment_genes_clint, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Engulfment_genes_clint, showlabels = F, Coi=c("")) + 
  ggtitle("Engulfment of apoptotic cell") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Engulfment_genes_clint, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Engulfment_genes_clint, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_Engulfment_genes_clint_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Engulfment_genes_clint, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_Engulfment_genes_clint_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(Engulfment_genes_clint, Gdiffexp)
nonSigGenes<-setdiff(Engulfment_genes_clint, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Engulfment_genes_clint_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Engulfment_genes_clint_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Engulfment_genes_clint_nonExp", append = T)
}


########################################
##### Phagocytosis, engulfment (Clint) -> Full, with children terms
########################################

#Read in Amigo results from excel file
Phagocytosis_engulfment <- read_xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/documentation/Full_GO_phagocytosis, engulfment_0006911.xlsx", col_names = F)

Phagocytosis_engulfment_genes<-unique(Phagocytosis_engulfment[[1]])


setdiff(Phagocytosis_engulfment_genes, rownames(expTable_mean))
setdiff(Phagocytosis_engulfment_genes, rownames(countData))

nonExpGenes <-setdiff(Phagocytosis_engulfment_genes, rownames(expTable_mean))
Phagocytosis_engulfment_genes<-intersect(Phagocytosis_engulfment_genes, rownames(expTable_mean))

p1<-plotDotplot(barycoords, Gdiffexp, Goi=Phagocytosis_engulfment_genes, showlabels = F, Coi=c("")) + 
  ggtitle("Phagocytosis/Engulfment genes") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2<-plotRoseplot(barycoords, Gdiffexp, Goi=Phagocytosis_engulfment_genes, showlabels = F, Coi=c(""))
grid.arrange(p1, p2, ncol=2)

### Export plots
p<-plotDotplot(barycoords, Gdiffexp, Goi=Phagocytosis_engulfment_genes, showlabels = F) + 
  theme(legend.position = "none")
ggsave(p,file="../results/triwiseResults/plots/triwisePlot_Phagocytosis_engulfment_genes_final_paper.png",dpi=300)

p<-plotRoseplot(barycoords, Gdiffexp, Goi=Phagocytosis_engulfment_genes, showlabels = F)
ggsave(p,file="../results/triwiseResults/plots/rosePlot_Phagocytosis_engulfment_genes_final_paper.png",dpi=300)


### Split in DE and non-DE genes
sigGenes<-intersect(Phagocytosis_engulfment_genes, Gdiffexp)
nonSigGenes<-setdiff(Phagocytosis_engulfment_genes, Gdiffexp)
matrixSigGenes<-barycoords[sigGenes,] %>% .[order(.$r, decreasing = T),]
matrixNonSigGenes<-barycoords[nonSigGenes,] %>% .[order(.$r, decreasing = T),]
write.xlsx(matrixSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_DE", append = T)
write.xlsx(matrixNonSigGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_nonDE", append = T)
if (length(nonExpGenes) > 1) {
  write.xlsx(nonExpGenes, "../results/heatmaps/Final_paper_Simon/summary_ownGOlists_final_paper_oct2019.xls", "Apoptotic_Clearance_genes_nonExp", append = T)
}


### Overlap lists
write.xlsx(Engulfment_genes_clint, "../results/heatmaps/Final_paper_Simon/overlap_engulfment_gene_lists_full.xls", "Engulfment_genes_clint")
write.xlsx(engulfmentGenes, "../results/heatmaps/Final_paper_Simon/overlap_engulfment_gene_lists_full.xls", "Engulfment_genes", append = T)
write.xlsx(Apop_cell_clearance_genes, "../results/heatmaps/Final_paper_Simon/overlap_engulfment_gene_lists_full.xls", "Apop_cell_clearance_genes", append = T)
write.xlsx(Phagocytosis_engulfment_genes, "../results/heatmaps/Final_paper_Simon/overlap_engulfment_gene_lists_full.xls", "Phagocytosis_engulfment_genes", append = T)

######################################################################
########## TEST SIGNIFICANCE OF OWN GO TERMS
######################################################################
gsetsNew<-gsets
gsetindexNew<-gsetindex


gsetsNew$`UPR`<-UPRgenesList2
gsetsNew$`Xbp1`<-Xbp1Genes_v1
gsetsNew$`IsrAtf`<-IsrAtf
gsetsNew$`atf6`<-atf6
gsetsNew$`riddHan`<-riddHan
gsetsNew$`riddGerlo`<-riddGerlo
gsetsNew$`riddMerge`<-riddMerge
gsetsNew$`lipidGenes`<-lipidGenes
gsetsNew$`cholesterolGenes`<-cholesterolGenes
gsetsNew$`Amigo_chol_genes_known`<-Amigo_chol_genes_known
gsetsNew$`effluxGenes`<-effluxGenes 
gsetsNew$`biosyntheseGenes`<-biosyntheseGenes 
gsetsNew$`cholesterolGenesAll`<-cholesterolGenesAll 
gsetsNew$`DCmaturationTLR`<-DCmaturationTLR
gsetsNew$`DCmaturationHomeo`<-DCmaturationHomeo
gsetsNew$`DCmaturationCommon`<-DCmaturationCommon
gsetsNew$`engulfmentGenes`<-engulfmentGenes 

#Adapted lines by Clint
gsetindexNew<-rbind(gsetindexNew,rep("UPR",4), rep("Xbp1",4), rep("IsrAtf",4), rep("atf6",4), rep("riddHan",4), rep("riddGerlo",4), rep("riddMerge",4), 
                    rep("lipidGenes",4), rep("cholesterolGenes",4), rep("Amigo_chol_genes_known",4), rep("effluxGenes",4), rep("biosyntheseGenes",4),
                    rep("cholesterolGenesAll",4), rep("DCmaturationTLR",4), rep("DCmaturationHomeo",4), rep("DCmaturationCommon",4), rep("engulfmentGenes",4))
wantedNames<-tail(gsetindexNew$name,17)


## Using rank statistic:
scoresTest = testUnidirectionality(barycoords, gsetsNew, Gdiffexp, statistic = "rank", nsamples=1e+6) #Not on windows: mc.cores = 8,
saveRDS(scoresTest,file="../results/scores_GOterms_test_new3.rds") #Version 1 (10), version 2 (11), version 3 (17)
# scoresTest <- readRDS("../results/scores_GOterms_test_new3.rds")

## Using diffexp statistic:
scoresTest = testUnidirectionality(barycoords, gsetsNew, Gdiffexp, statistic = "r", nsamples=1e+6) #Not on windows: mc.cores = 8,
saveRDS(scoresTest,file="../results/scores_GOterms_test_r.rds")
# scoresTest <- readRDS("../results/scores_GOterms_test_r.rds")

## Using diffexp statistic:
scoresTest = testUnidirectionality(barycoords, gsetsNew, Gdiffexp, statistic = "diffexp", nsamples=1e+6) #Not on windows: mc.cores = 8,
saveRDS(scoresTest,file="../results/scores_GOterms_test_diffexp.rds")
# scoresTest <- readRDS("../results/scores_GOterms_test_diffexp.rds")


scoresTestRes = left_join(scoresTest, gsetindexNew, by="gsetid") %>% filter(qval < 0.05) %>% arrange(qval, z)

scoresTestRes = scoresTestRes[(scoresTestRes$qval < 0.05) & (scoresTestRes$z > 3500), ] #rank
scoresTestRes = scoresTestRes[(scoresTestRes$qval < 0.05) & (scoresTestRes$z > 0.15), ] #r statistic
scoresTestRes = scoresTestRes[(scoresTestRes$qval < 0.05) & (scoresTestRes$z > 0.05), ] #diffexp

dim(scoresTestRes)
##107 with extra cutoff z>3500!! 298 with cutoff z>2000!! (rank)
##129 with cutoff z>0.15!! (r)
##135 with cutoff z>0.05!! (diffexp)

### List all terms
tail(scoresTest,17)
### List significant terms
scoresTestRes[scoresTestRes$name %in% wantedNames,]

write.xlsx(scoresTestRes[scoresTestRes$name %in% wantedNames,], file = "../results/triwiseResults/Significance_GO_new_diffexp.xlsx")

######################################################################
########## CREATE HEATMAP
######################################################################

matrixAllDE<-rbind(DEgenesGroup1,DEgenesGroup2,DEgenesGroup3)
matrixAllDE <-matrixAllDE[unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2),rownames(DEgenesGroup3))),]

##Adjustments from Sophie:
adjustedCholesterol<-c(intersect(cholesterolGenes, Gdiffexp),"Nus1","Lss","Tm7sf2","Sc4mol","Nsdhl","Dhcr24","Sc5d","Dhcr7","Ebp") #Include non-DE genes
adjustedDCmaturationHomeo<-rownames(matrixAllDE[intersect(DCmaturationHomeo, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
adjustedDCmaturationHomeo<-c(adjustedDCmaturationHomeo,c("Nuak1","Pik3r3","Hecw2","Kcnn1","Scube3","Cdcp1","Ccdc120")) #Include non-DE genes
adjustedAmiGOChol <- setdiff(Amigo_chol_genes_known, c("Ccl3", "Nfkbia","Xbp1")) #Leave out. Confusing.
adjustedXBP1 <- c(Xbp1Genes_v1,c("Txndc5", "Dnajc10", "Sec61a1")) 

## Extra check adjusted genes: (Clint)
genes_chol_adjustment <- c("Nus1","Lss","Tm7sf2","Sc4mol","Nsdhl","Dhcr24","Sc5d","Dhcr7","Ebp")
datalistadjustment <- list()
datalistadjustment2 <- list()
datalistadjustment3 <- list()

for (i in 1:length(genes_chol_adjustment)) {
  if (genes_chol_adjustment[i] %in% rownames(allGenesGroup1)){
    datalistadjustment[[i]] <- allGenesGroup1[genes_chol_adjustment[i],]
  }
  if (genes_chol_adjustment[i] %in% rownames(allGenesGroup2)){
    datalistadjustment2[[i]] <- allGenesGroup2[genes_chol_adjustment[i],]
  }
  if (genes_chol_adjustment[i] %in% rownames(allGenesGroup3)){
    datalistadjustment3[[i]] <- allGenesGroup3[genes_chol_adjustment[i],]
  }
}

big_data = do.call(rbind, datalistadjustment) %>% .[order(.$adj.P.Val, decreasing = F),]
big_data2 = do.call(rbind, datalistadjustment2) %>% .[order(.$adj.P.Val, decreasing = F),]
big_data3 = do.call(rbind, datalistadjustment3) %>% .[order(.$adj.P.Val, decreasing = F),]

write.xlsx(big_data, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_chol_1")
write.xlsx(big_data2, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_chol_2", append = T)
write.xlsx(big_data3, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_chol_3", append = T)

genes_DChomeo_adjustment <- c("Nuak1","Pik3r3","Hecw2","Kcnn1","Scube3","Cdcp1","Ccdc120")
datalistadjustment <- list()
datalistadjustment2 <- list()
datalistadjustment3 <- list()

for (i in 1:length(genes_DChomeo_adjustment)) {
  if (genes_DChomeo_adjustment[i] %in% rownames(allGenesGroup1)){
    datalistadjustment[[i]] <- allGenesGroup1[genes_DChomeo_adjustment[i],]
  }
  if (genes_DChomeo_adjustment[i] %in% rownames(allGenesGroup2)){
    datalistadjustment2[[i]] <- allGenesGroup2[genes_DChomeo_adjustment[i],]
  }
  if (genes_DChomeo_adjustment[i] %in% rownames(allGenesGroup3)){
    datalistadjustment3[[i]] <- allGenesGroup3[genes_DChomeo_adjustment[i],]
  }
}

big_data = do.call(rbind, datalistadjustment) %>% .[order(.$adj.P.Val, decreasing = F),]
big_data2 = do.call(rbind, datalistadjustment2) %>% .[order(.$adj.P.Val, decreasing = F),]
big_data3 = do.call(rbind, datalistadjustment3) %>% .[order(.$adj.P.Val, decreasing = F),]

write.xlsx(big_data, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_DChomeo_1", append = T)
write.xlsx(big_data2, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_DChomeo_2", append = T)
write.xlsx(big_data3, "../results/triwiseResults/Expression_adjusted_genelists.xls", "adjusted_DChomeo_3", append = T)

###Get genelistname
genelistname <- "adjustedAmiGOcholesterolv3"
genelistname <- "adjustedXBP1"
genelistname <- "DCmaturationHomeo"
genelistname <- "DCmaturationCommon"
genelistname <- "RIDDmerge"
genelistname <- "Isr"


### Get genes
wantedGenes<-rownames(matrixAllDE[intersect(adjustedAmiGOChol, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(adjustedXBP1, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(adjustedDCmaturationHomeo, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(DCmaturationCommon, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(DCmaturationCommon, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),]) %>% tail(.,30)
wantedGenes<-rownames(matrixAllDE[intersect(riddMerge, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(Isr, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),])
wantedGenes<-rownames(matrixAllDE[intersect(Isr, Gdiffexp),] %>% .[order(.$adj.P.Val, decreasing = T),]) %>% tail(.,30)

### Get samples
colsWT<-c(grep("WT",colnames(expTable)), grep("fl_fl",colnames(expTable)))
colsIre1KO<-grep("IRE1_KO",colnames(expTable))
colsXpb1KO<-grep("XBP1_KO",colnames(expTable))

wantedSamples<-c(colsWT,colsXpb1KO,colsIre1KO)

### Normalize
expProfiles<-expTable[wantedGenes,wantedSamples]

expProfilesNorm<-normalizePerGene(expProfiles)
gem<-apply(expProfiles,1,mean)
expProfilesNorm2<-expProfiles-gem

# Get order according to angle (starts at -3 and goes to 3: counterclockwise!!)
expProfilesNorm_Angle_ordered <- barycoords[rownames(expProfilesNorm2),]
expProfilesNorm_Angle_ordered <- expProfilesNorm_Angle_ordered[order(expProfilesNorm_Angle_ordered[,3]),]
expProfilesNorm2<-expProfilesNorm2[rownames(expProfilesNorm_Angle_ordered),]
#expProfilesNorm2 <- cbind(expProfilesNorm2, barycoords[rownames(expProfilesNorm2),"angle"])


## Change column names heatmap
colnames(expProfilesNorm2) <- c(rep("XBP1/IRE1 fl/fl",3), rep(paste0("XBP1 ", as.character("\u0394")),3), rep(paste0("XBP1/IRE1 ", as.character("\u0394")),3))

### Prepare heatmap and Triwise

##Create colors
myColorPalette<-c("#08306b", "#08326e", "#083573", "#083876", "#083b7b", "#083d7f", "#084083", "#084387", "#08468c", "#084990", 
                  "#084b94", "#0a4f97", "#0b5199", "#0e559d", "#0f579f", "#125ba2", "#135da4", "#1661a7", "#1864aa", "#1967ad", 
                  "#1c6ab0", "#1f6db1", "#2372b4", "#2675b6", "#2a79b8", "#2e7ebb", "#3181bc", "#3685bf", "#3888c1", "#3d8dc3", 
                  "#4090c5", "#4794c7", "#4c98c9", "#539dcb", "#5ba1ce", "#60a5d0", "#67aad2", "#6cadd4", "#74b2d6", "#79b5d8", 
                  "#81badb", "#8bbfde", "#99c7e2", "#a8d0e6", "#b2d5e9", "#c1dded", "#cbe3f0", "#daebf4", "#e4f0f7", "#f3f8fb", 
                  "#fffefd", "#fff8f5", "#feefe9", "#fee9e1", "#fee0d4", "#fedacc", "#fdd1c0", "#fdcbb8", "#fdc2ac", "#fcb9a0", 
                  "#fcb398", "#fcab8f", "#fca68a", "#fc9e82", "#fc997c", "#fc9174", "#fc8c6e", "#fb8466", "#fb7d5d", "#fb7758", 
                  "#fb7050", "#f96a4d", "#f66348", "#f45e45", "#f15640", "#ee4e3b", "#ec4937", "#ea4133", "#e83c2f", "#e5342a", 
                  "#e32f27", "#dd2c25", "#d92924", "#d32622", "#cd2220", "#c9201f", "#c31d1d", "#bf1a1c", "#b9171a", "#b51419", 
                  "#ae1117", "#a91016", "#a00e15", "#980c14", "#920a13", "#890812", "#840711", "#7b0510", "#75030f", "#6d010e")

paletteLength<-100

## Extra code annotation rows
Origin_genes<-c(rep("Original list",6),"Manually added","Original list","Manually added",rep("Original list",18), rep("Manually added",2),rep("Original list",2))
df_genes2 <- data.frame(wantedGenes,Origin_genes)
rownames(df_genes2) <- df_genes2$wantedGenes
df_genes2 <- df_genes2[,-1, drop=F]
my_color2 <- list(Origin_genes = c("Original list" = "red", "Manually added" = "blue"))
                              
### Create heatmap
myBreaks <- c(seq(min(expProfilesNorm2), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(expProfilesNorm2)/paletteLength, max(expProfilesNorm2), length.out=floor(paletteLength/2)))

p<-pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 11, cellheight = 11, gaps_col = c(3,6), gaps_row = c(21,24,29), fontsize=12,
         show_rownames = T, show_colnames = T, breaks = myBreaks)

# p<-pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
#             border_color = "gray25", cellwidth = 11, cellheight = 11, gaps_col = c(3,6), gaps_row = c(3,6), fontsize=12,
#             show_rownames = T, show_colnames = T, breaks = myBreaks, annotation_row = df_genes2, annotation_colors = my_color2)

print(p)

ggsave(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_heatmap.png"),dpi=300, height = 10)


cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_heatmap.pdf"),height=10)
pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 11, cellheight = 11, gaps_col = c(3,6), gaps_row = c(21,24,29), fontsize=12,
         show_rownames = T, show_colnames = T, breaks = myBreaks)
# pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
#          border_color = "gray25", cellwidth = 11, cellheight = 11, gaps_col = c(3,6), gaps_row = c(3,6), fontsize=12,
#          show_rownames = T, show_colnames = T, breaks = myBreaks, annotation_row = df_genes2, annotation_colors = my_color2)

dev.off()

### Create triwise
##Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE, Goi = wantedGenes) #top20 indicated with line. Gpin doesn't work and color doesn't either
print(p)
saveWidget(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_triwisePlot_withLines.html")) ##needs full path!

###Triwise plot (red color) (choose correct genelist!!)
p<-plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F) + 
  theme(legend.position = "none") #Can't use wantedGenes!! Must use list before intersect with diffexp!!
print(p)
ggsave(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_triwisePlot_red.png"),dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_triwisePlot_red.pdf"))
plotDotplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F) + 
  theme(legend.position = "none") #Can't use wantedGenes!! Must use list before intersect with diffexp!!
dev.off()

##Rose plot (choose correct genelist)
p<-plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F) ##Use genelist before intersect (but not necessary)!
print(p)
ggsave(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_rosePlot.png"),dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/heatmaps/Final_paper_Simon/",genelistname,"_rosePlot.pdf"))
plotRoseplot(barycoords, Gdiffexp, Goi=DCmaturationCommonNew, showlabels = F) ##Use genelist before intersect (but not necessary)!
dev.off()
