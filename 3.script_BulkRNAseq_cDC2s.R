## Script for analyzing Xbp1/Ire1KO cDC2 bulk RNA-seq project data
## limma-edgeR DEG pipeline 

library("limma")
library("edgeR")
library("ggplot2")
library("triwise")
library("htmlwidgets")
library("rgl")
library("RColorBrewer")
library("gplots")
library('xlsx')

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

### Only cDC2 data without XBP1_IRE1_KO_cDC2_3 and XBP1_KO_cDC2_4 (due to low lib size)
colData<-colData[c(13:18,20:23),]
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
theColors<-c(rep("royalblue",4),rep("darkgreen",3),rep("limegreen",3))

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

###########################################
###########################################

## Update 2023 for upload to GEO

# Check names
all(colnames(countData) == colnames(expTable))

write.table(countData, "../results/GEO/Bulk_RNA_seq_raw_counts_cDC2s_Ire1_Xbp1_experiment.txt", sep="\t")
write.table(expTable, "../results/GEO/Bulk_RNA_seq_normalized_counts_cDC2s_Ire1_Xbp1_experiment.txt", sep="\t")

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

hmcol <- colorRampPalette(brewer.pal(10, "GnBu"))(100)

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
cont.matrix <- makeContrasts(group1=cDC2.XBPKO-cDC2.WT, group2=cDC2.IREKO-cDC2.WT, group3=cDC2.IREKO-cDC2.XBPKO, levels=design)
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
##10

##############################
##### Ire1 KO vs WT
##############################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,1)
dim(DEgenesGroup2)
##3

##############################
##### Ire1 KO vs Xbp1 KO
##############################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.05,1)
dim(DEgenesGroup3)
##3

### Write results
write.xlsx(DEgenesGroup1, "../cDC2_analysis_Clint/summary_DEgenes_cDC2.xls", "Xbp1KO_vs_WT")
write.xlsx(DEgenesGroup2, "../cDC2_analysis_Clint/summary_DEgenes_cDC2.xls", "Ire1KO_vs_WT", append=T)
write.xlsx(DEgenesGroup3, "../cDC2_analysis_Clint/summary_DEgenes_cDC2.xls", "Ire1KO_vs_Xbp1KO", append=T)

####################################################################################################
############################## ----------- TRIWISE PLOTS ------------ ##############################
####################################################################################################

##### All DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2),rownames(DEgenesGroup3)))
length(allDEgenes)
##10

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

saveRDS(barycoords, "../results/barycoords_cDC2s.rds")
saveRDS(Gdiffexp, "../results/Gdiffexp_cDC2s.rds")

###Triwise plot
p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors7, showlabels = F, sizevalues = 50)
print(p)
ggsave(p,file="../cDC2_analysis_Clint/triwisePlot.png",dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/cDC2_analysis_Clint/triwisePlot.pdf"))
plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors7, showlabels = F, sizevalues = 50)
dev.off()

###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)
ggsave(p,file="../cDC2_analysis_Clint/rosePlot.png",dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/cDC2_analysis_Clint/rosePlot.pdf"))
plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F) 
dev.off()

###Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE, Goi = allDEgenes)
print(p)
saveWidget(p,file="/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/cDC2_analysis_Clint/interactiveTriwisePlot.html") ##needs full path!

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
ggsave(p,file="../cDC2_analysis_Clint/triwisePlot_withColors_newcolors.png",dpi=300)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/cDC2_analysis_Clint/triwisePlot_withColors_newcolors.pdf"))
plotDotplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F, Goi=list(Set1=genesSet1,Set2=genesSet2,Set3=genesSet3,Set4=genesSet4,
                                                                         Set5=genesSet5,Set6=genesSet6,Set7=genesSet7), colorvalues=wantedColors7) +
  theme(legend.position = "none")
dev.off()

myList<-list("set1"=genesSet1,"set2"=genesSet2,"set3"=genesSet3,"set4"=genesSet4,"set5"=genesSet5,"set6"=genesSet6,"set7"=genesSet7)
saveRDS(myList,file="../cDC2_analysis_Clint/listGenes_coloredTriwise_cDC2s.rds")
