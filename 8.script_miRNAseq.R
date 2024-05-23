## Script for analyzing Xbp1/Ire1KO cDC1 miRNA-seq project data
## limma-edgeR DEmiRNA pipeline combined with triwise analysis

library("limma")
library("edgeR")
library("ggplot2")
library('miRBaseConverter')
library('tidyverse')
library("rgl")
library("RColorBrewer")
library("gplots")
library('openxlsx')
library("triwise")
library("htmlwidgets")
library("org.Mm.eg.db")

##### Xbp1_WT + Ire1Xbp1_WT = 8 replicates
##### Xbp1_KO = 4 replicates
##### Ire1Xbp1_KO = 4 replicates

setwd("/home/clintdn/VIB/DATA/Sophie/miRNASeq/")

dir.create("results_Clint_2023")
dir.create("results_Clint_2023/Robjects")
dir.create("results_Clint_2023/triwiseResults")

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

###Get DE genes, working on uncorrected P-val
getDEgenes_uncorr<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$P.Value<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}

## More functions!!

# barycoords=theBarycoords
# Gdiffexp = allDEgenes_cDC2
# Goi = NULL 
# Coi = attr(barycoords, "conditions")
# colorby = "diffexp" 
# colorvalues = wantedColors
# rmax = 5
# showlabels = TRUE
# sizevalues = stats::setNames(c(0.5,2), c(FALSE, TRUE))
# alphavalues = stats::setNames(c(0.8,0.8), c(FALSE, TRUE))
# barycoords2 = NULL
# baseangle = 0

plotDotplot_byLM<-function (barycoords, Gdiffexp = rownames(barycoords), Goi = NULL, 
                            Coi = attr(barycoords, "conditions"), colorby = "diffexp", 
                            colorvalues = NULL, rmax = 5, showlabels = TRUE, sizevalues = stats::setNames(c(0.5, 
                                                                                                            2), c(FALSE, TRUE)), alphavalues = stats::setNames(c(0.8, 
                                                                                                                                                                 0.8), c(FALSE, TRUE)), barycoords2 = NULL, baseangle = 0) 
{
  if (!is.list(Goi)) {
    Goi = list(gset = Goi)
  }
  barypoints = as.data.frame(barycoords)
  barypoints$diffexp = rownames(barypoints) %in% Gdiffexp
  barypoints$ingset = FALSE
  barypoints$gsetname = "all"
  if (is.null(names(Goi))) 
    names(Goi) = 1:length(Goi)
  for (gsetname in names(Goi)) {
    barypoints[Goi[[gsetname]], "ingset"] = TRUE
    barypoints[Goi[[gsetname]], "gsetname"] = gsetname
  }
  barypoints$type = paste0(ifelse(barypoints$diffexp, "diff", 
                                  "nodiff"), barypoints$gsetname)
  if (colorby == "diffexp") {
    if (is.null(colorvalues)) {
      colorvalues = stats::setNames(c("#333333", RColorBrewer::brewer.pal(max(length(Goi), 
                                                                              3), "Set1")), c("diffall", paste0("diff", names(Goi))))
      colorvalues = c(colorvalues, stats::setNames(c("#AAAAAA", 
                                                     RColorBrewer::brewer.pal(max(length(Goi), 3), 
                                                                              "Pastel1")), c("nodiffall", paste0("nodiff", 
                                                                                                                 names(Goi)))))
    }
    color = ggplot2::scale_colour_manual(values = colorvalues, 
                                         name = "type")
    barypoints$colorby = barypoints$type
    barypoints$alphaby = barypoints$diffexp
    barypoints$sizeby = barypoints$diffexp
  } else if (colorby == "z") {
    if (!("z" %in% colnames(barypoints))) 
      stop("z column not defined")
    color = ggplot2::scale_colour_continuous()
    barypoints$colorby = barypoints$z
    barypoints$alphaby = TRUE
    barypoints$sizeby = TRUE
  }  else {
    color = ggplot2::scale_colour_continuous()
    barypoints$colorby = barypoints$r
    barypoints$alphaby = TRUE
    barypoints$sizeby = TRUE
  }
  alpha = ggplot2::scale_alpha_manual(values = stats::setNames(c(0.8, 
                                                                 0.8), c(FALSE, TRUE)), name = "diffexp")
  size = ggplot2::scale_size_manual(values = stats::setNames(c(0.8, 
                                                               1.5), c(FALSE, TRUE)), name = "diffexp")
  order = with(barypoints, order(ingset, diffexp))
  plot = drawHexagonGrid(rmax, showlabels = showlabels, baseangle = baseangle) + 
    drawDirections(rmax, Coi, baseangle = baseangle) + drawDotplot(barypoints, 
                                                                   rmax, color = color, order = order, alpha = alpha, size = size, 
                                                                   baseangle = baseangle)
  if (!is.null(barycoords2) && sum(barypoints$ingset) > 0) {
    barypoints2 = as.data.frame(barycoords2)[rownames(barypoints), 
                                             ]
    plot = plot + drawConnectionplot(barypoints[barypoints$ingset, 
                                                ], barypoints2[barypoints$ingset, ], rmax, baseangle = baseangle)
  }
  plot
}


drawHexagonGrid <- function(rmax=5, color="#999999", showlabels=TRUE, baseangle=0) {
  radii = seq_len(rmax)
  hexagonpoints <- plyr::ldply(seq(0, pi*2-0.001, pi/3), function(angle) data.frame(x=cos(angle+baseangle), y=sin(angle+baseangle)))
  gridData <- dplyr::bind_rows(lapply(radii, function(r) data.frame(r=r, x=hexagonpoints$x * r, y =hexagonpoints$y * r)))
  
  plot =  drawGridBasis() +
    ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, group=r), gridData, fill=NA, colour=color) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1))
  
  if (showlabels) {
    labelData <- data.frame(r=radii, label=2^radii, y=sapply(radii, function(x) hexagonPolar(pi/2, x)), x=0)
    plot = plot + ggplot2::geom_label(ggplot2::aes(x=x, y=y, label=label), hjust=0.5, data=labelData)
  }
  plot
}

drawGridBasis <- function() {
  ggplot2::ggplot() +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x=ggplot2::element_blank(),
      axis.title.y=ggplot2::element_blank(),
      axis.ticks=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_blank(),
      panel.background=ggplot2::element_blank()
    )
}

drawDirections <- function(rmax, labels, labelmargin=0.05, type="hexagonal", baseangle=0) {
  directionsData = plyr::ldply(seq(0, pi*2-0.001, pi/3*2), function(angle) {
    angle = angle + baseangle
    x = cos(angle) * rmax
    y = sin(angle) * rmax
    
    if (y<0){
      va=1
    } else if (y == 0) {
      va=0.5
    } else{
      va=0
    }
    if (x<0){
      ha = 1
    } else if (x < 0.1 & x > -0.1) {
      ha=0.5
    } else{
      ha=0
    }
    
    rot = 0
    if(sin(angle) >= 0) {
      rot = angle - pi/2
      va = 0
      ha = 0.5
    } else {
      rot = angle + pi/2
      va = 1
      ha = 0.5
    }
    
    rot = rot / pi * 180
    
    data.frame(x=x, y=y, va=va, ha=ha, xlabel=x*(1+labelmargin),ylabel=y*(1+labelmargin), rot=rot)
  })
  directionsData$label = labels
  
  c(
    ggplot2::geom_text(ggplot2::aes(label=label, x=xlabel, y=ylabel, vjust=va, hjust=ha, angle=rot), data=directionsData),
    ggplot2::geom_segment(ggplot2::aes(xend=0, x=x, yend=0, y=y, group=label), data=directionsData, alpha=0.5)
  )
}

drawDotplot <- function(barypoints, rmax=5, color=ggplot2::scale_color_grey(), alpha=ggplot2::scale_alpha(), size=ggplot2::scale_size(), order=NULL, baseangle=0) {
  barypoints = clipHexagon(barypoints, rmax, baseangle)
  
  if(!is.null(order)) {
    barypoints = barypoints[order,]
  }
  
  dotplot = ggplot2::geom_point(ggplot2::aes(x=xclip, y=yclip, color=colorby, alpha=alphaby, size=sizeby), data=barypoints)
  dotplot = c(dotplot, color, alpha, size)
  
  dotplot
}

drawConnectionplot <- function(barypoints, barypoints2, rmax=5, order=NULL, baseangle=0) {
  barypoints = clipHexagon(barypoints, rmax, baseangle)
  barypoints2 = clipHexagon(barypoints2, rmax, baseangle)
  
  colnames(barypoints2) = paste0(colnames(barypoints2), "2")
  
  allbarypoints = data.frame(barypoints, barypoints2)
  
  if(!is.null(order)) {
    allbarypoints = allbarypoints[order,]
  }
  
  dotplot = ggplot2::geom_segment(ggplot2::aes(x=xclip, y=yclip, xend=xclip2, yend=yclip2, color=colorby), data=allbarypoints, arrow=ggplot2::arrow(type="closed", length=ggplot2::unit(0.05, "inches")))
  
  dotplot
}

clipHexagon <- function(barypoints, rmax, baseangle=0) {
  barypoints["rclip"] = mapply(function(angle, r) {min(hexagonPolar(angle, rmax, baseangle), r)}, barypoints$angle, barypoints$r)
  barypoints["xclip"] = cos(barypoints["angle"]+baseangle) * barypoints["rclip"]
  barypoints["yclip"] = sin(barypoints["angle"]+baseangle) * barypoints["rclip"]
  
  barypoints
}


hexagonPolar <- function(angle, radius=1, baseangle=0) {
  delta <- 2*pi/6
  cos(delta/2)/cos((angle %% delta)-delta/2) * radius
}

################################################################################
######### LOAD DATA
################################################################################

getwd()

##### Step1: get all file names #####
dirNames <- list.dirs(path = "./ruwe data Pieter Mestdagh/rerun analysis/rerun_2510/", full.names=FALSE)
##Remove first directory = data directory
dirNames<-dirNames[-1]
countFiles<-paste0(dirNames,"/",dirNames,"_miRs.txt")

##### Step2: parse all the files #####
DT <- list()
all_miRNA<-c()

for(fileName in countFiles){
  countData<-read.table(file=paste0("ruwe data Pieter Mestdagh/rerun analysis/rerun_2510/",fileName),sep="\t", header=FALSE, stringsAsFactors=FALSE)
  countData<-countData[,c(2,3,1)]
  countData<-countData[,-3]
  
  tmp<-strsplit(fileName,"/")[[1]][1]
  colnames(countData)<-c("ID",tmp)
  
  DT[[fileName]]<-countData
  all_miRNA<-c(all_miRNA, countData$ID)
}

all_miRNA<-all_miRNA[! duplicated(all_miRNA)]
length(all_miRNA)

##### Step3: make sure that everyone has all miRNAs #####
newDT<-lapply(DT, function(x){
  missing_miRNA<-setdiff(all_miRNA, x$ID)
  toAdd<-data.frame(missing_miRNA,missing_miRNA)
  colnames(toAdd)<-colnames(x)
  toAdd[,2]<-0
  new<-rbind(x, toAdd)
  x<-new
})

##### Step4: merge all files #####
countData <- newDT[[countFiles[1]]]

##Add the other count files
for (i in 2:length(countFiles)) {
  y <- newDT[[countFiles[i]]]
  z <- merge(countData, y, by = c("ID"))
  countData <- z
}

##### Clean up table #####
rownames(countData)<-countData[,1]
countData<-countData[,-1]
dim(countData)
# 853  20

##### Load meta data #####
colData<-read.table(file="Liesbet_2021/metadata.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE)
rownames(colData)<-colData$fileName
colData<-colData[,-1]
dim(colData)
# 20  2

##### Change colnames #####
converter<-read.table(file="documentation/convert_sampleNames.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE)
countData<-countData[,converter$oldName]
cbind(colnames(countData),converter)
colnames(countData)<-converter$newName

dim(countData)
# 853  20

##### Reorder columns #####
cols_WT<-grep('WT_',colnames(countData)) #Change in analysis
cols_Xbp1_KO<-grep('^Xbp1_KO_',colnames(countData))
cols_Ire1_KO<-grep('^Ire1Xbp1_KO_',colnames(countData))

countData<-countData[,c(cols_WT,cols_Xbp1_KO,cols_Ire1_KO)]
dim(countData)
# 853  20

## Fix coldata -> same order!!
colData<-colData[colnames(countData),]
dim(colData)
cbind(colnames(countData),rownames(colData))

## Fix coldata: new condition column (combine Xbp1WT and Ire1WT)
colData$condition_new<-as.character(colData$condition)
colData$condition_new[1:12]<-"WT"
colData$condition_new<-as.character(colData$condition_new)

## Remove cDC2 samples
##### Remove Xbp1_WT_1 (5_mouse6) #####
theID<-which(colnames(countData)=="Xbp1_WT_1")
countData<-countData[,-theID]
colData<-colData[-theID,]
dim(countData)
# 853  19

##### Remove 7_mouse7, 11_mouse4, 19_mouse3) #####
theIDs<-c(which(colnames(countData)=="Xbp1_WT_2"),
          which(colnames(countData)=="Xbp1_WT_3"),
          which(colnames(countData)=="Xbp1_WT_7"))

countData<-countData[,-theIDs]
colData<-colData[-theIDs,]
dim(countData)
# 853  16

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

keep = rowSums(cpm(y)>1) >= 4
y = y[keep,]
dim(y)
##473  16
dim(yNoFilter)
##853  16

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
##Remove cDC2
theColors<-c(rep("royalblue",4),rep("darkgreen",4),rep("darkorange",4),rep("darkorchid",4))
cbind(colData, theColors)


png(file="results_Clint_2023/4_MDSplot.png", width = 1515, height = 1138)
plotMDS(y, cex=1.2, col=theColors)
dev.off()

################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
TS <- paste(colData$cell, colData$condition_new, sep=".")
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
png(file="results_Clint_2023/meanVariancePlot.png",width = 1515, height = 1138)
v <- voom(y, design, plot = TRUE)
dev.off()

expTable<-v$E

saveRDS(expTable, "results_Clint_2023/Robjects/expTable_2023.rds")

##### Normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

################################################################################
########## CHECK FILTERING AND NORMALISATION
################################################################################

#################### BARPLOT ####################

##### Barplot lib sizes raw counts
png(file="results_Clint_2023/1_barplot_beforeNorm.png", width = 1515, height = 1138)
par(mar = c(9,3,3,1)) #more margin: bottom, left, top, right
bp<-barplot(yNoFilter$samples$lib.size*1e-6,axisnames=FALSE,main="Barplot lib sizes of raw counts",ylab="Library size (millions)")
axis(1, labels=rownames(yNoFilter$samples), at = bp, las=2, cex.axis=0.8)
dev.off()

##### Barplot lib sizes normalised counts
# countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

png(file="results_Clint_2023/1_barplot_afterNorm.png", width = 1515, height = 1138)
par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
bp<-barplot(colSums(countData_norm)*1e-6,axisnames=FALSE,main="Barplot lib sizes of normalised counts",ylab="Library size (millions)")
axis(1, labels=colnames(countData_norm), at = bp, las=2, cex.axis=0.7)
dev.off()

#################### BOXPLOT ####################
col <- rainbow(nrow(colData))

y2<-y
y2$samples$norm.factors<-1
y2$counts[,1]<-ceiling(y2$counts[,1]*0.05)
y2$counts[,2]<-y2$counts[,2]*5

# par(mfrow=c(1,2), mar = c(12,4,3,1)) #more margin: bottom, left, top, right
png(file="results_Clint_2023/2_boxplot_beforeNorm.png", width = 1515, height = 1138)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Unnormalised data", ylab="log-cpm", col=col)
dev.off()

png(file="results_Clint_2023/2_boxplot_afterNorm.png", width = 1515, height = 1138)
y2<-calcNormFactors(y2)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Normalised data", ylab="log-cpm", col=col)
dev.off()

#################### DENSITY PLOT ####################
col <- topo.colors(nrow(colData))

png(file="results_Clint_2023/3_densityPlot.png", width = 1515, height = 1138)
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
dev.off()

par(mfrow=c(1,1))

#################### HISTOGRAM OF EXPTABLE ####################

png(file="results_Clint_2023/4_histogramFiltering.png", width = 1515, height = 1138)
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
dev.off()


par(mfrow=c(1,1))

################################################################################
########## PCA
################################################################################

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.99*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

PoV <- pca$sdev^2/sum(pca$sdev^2)
summary(pca)

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.5)


### Save 3D
dirToSave<-paste(getwd(), "/results_Clint_2023/",sep="")
writeWebGL(dir = dirToSave, filename = file.path(dirToSave, "4_pca.html"),
           template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
           snapshot = TRUE, font = "Arial")

### Save as image
rgl.viewpoint(0, 0)
rgl.snapshot("results_Clint_2023/4_pca_view1.png")
rgl.viewpoint(35, 0)
rgl.snapshot("results_Clint_2023/4_pca_view2.png")


################################################################################
########## CORRELATION HEATMAP SAMPLES
################################################################################

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

##heatmap 1: based on distance
distsRL <- dist(t(expTable),method="euclidean")
hc <- hclust(distsRL,method="ward.D")

pdf("results_Clint_2023/4_corrSamples_distance.pdf")
heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13), cexRow=0.9, cexCol=0.9)
dev.off()

##heatmap 2: based on correlation
cm=cor(expTable)

distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")

pdf("results_Clint_2023/4_corrSamples_correlation.pdf")
heatmap.2(cm, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=hmcol, margin=c(13, 13), cexRow=0.9, cexCol=0.9)
dev.off()


################################################################################
########## GET DE GENES
################################################################################

#### Fit linear model on data
fit <- lmFit(v, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=cDC1.Xbp1KO-cDC1.WT, 
                             group2=cDC1.Ire1KO-cDC1.WT, 
                             group3=cDC1.Ire1KO-cDC1.Xbp1KO,
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Quick list of DE genes
summa.fit <- decideTests(fit.eb, p.value = 0.01, lfc = 1)
summary(summa.fit)
# group1 group2 group3
# Down        2      1      0
# NotSig    470    470    473
# Up          1      2      0

summa.fit <- decideTests(fit.eb, p.value = 0.05, lfc = 1)
summary(summa.fit)
# group1 group2 group3
# Down        5      2      0
# NotSig    466    467    473
# Up          2      4      0

########################################
##### 1. Xbp1 KO vs WT
########################################
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.05,1)
DEgenesGroup1$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup1))[,2]
dim(DEgenesGroup1)
##7

DEgenesGroup1_clint<-getDEgenes(allGenesGroup1,0.15,0.5)
DEgenesGroup1_clint$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup1_clint))[,2]
dim(DEgenesGroup1_clint)
##25

DEgenesGroup1_lowS<-getDEgenes_uncorr(allGenesGroup1,0.05,1)
DEgenesGroup1_lowS$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup1_lowS))[,2]
dim(DEgenesGroup1_lowS)
##24

intersect(DEgenesGroup1_clint$miRNAname,DEgenesGroup1_lowS$miRNAname)

########################################
##### 2. Xbp1Ire1 KO vs WT
########################################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,1)
DEgenesGroup2$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup2))[,2]
dim(DEgenesGroup2)
##6

DEgenesGroup2_clint<-getDEgenes(allGenesGroup2,0.15,0.5)
DEgenesGroup2_clint$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup2_clint))[,2]
dim(DEgenesGroup2_clint)
##35

DEgenesGroup2_lowS<-getDEgenes_uncorr(allGenesGroup2,0.05,1)
DEgenesGroup2_lowS$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup2_lowS))[,2]
dim(DEgenesGroup2_lowS)
##28

########################################
##### 3. Xbp1Ire1 KO vs Xbp1 KO
########################################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.05,1)
DEgenesGroup3$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup3))[,2]
dim(DEgenesGroup3)
##0

DEgenesGroup3_clint<-getDEgenes(allGenesGroup3,0.15,0.5)
DEgenesGroup3_clint$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup3_clint))[,2]
dim(DEgenesGroup3_clint)
##18

DEgenesGroup3_lowS<-getDEgenes_uncorr(allGenesGroup3,0.05,1)
DEgenesGroup3_lowS$miRNAname = miRNA_AccessionToName(rownames(DEgenesGroup3_lowS))[,2]
dim(DEgenesGroup3_lowS)
##28

##### Put DE genes in list
coef<-c(1:ncol(cont.matrix))
listallgenes<-list()
listDEgenes<-list()
listDEgenes_lessStrict<-list()
listDEgenes_moreStrict<-list()

for(i in 1:length(coef)){
  allGenes<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=coef[i])
  DEgenes<-getDEgenes(allGenes,0.05,1)
  DEgenes$miRNAname = miRNA_AccessionToName(rownames(DEgenes))[,2]
  DEgenes_less_strict<-getDEgenes(allGenes,0.15,0.5)
  DEgenes_less_strict$miRNAname = miRNA_AccessionToName(rownames(DEgenes_less_strict))[,2]
  DEgenes_strict<-DEgenes[DEgenes$AveExpr >0,]
  allGenes$miRNAname = miRNA_AccessionToName(rownames(allGenes))[,2] #2022 added
  listallgenes[[i]]<-allGenes
  listDEgenes[[i]]<-DEgenes
  listDEgenes_lessStrict[[i]]<-DEgenes_less_strict
  listDEgenes_moreStrict[[i]]<-DEgenes_strict
}

names(listallgenes)<-c("cDC1.Xbp1KO-cDC1.WT", "cDC1.Ire1KO-cDC1.WT","cDC1.Ire1KO-cDC1.Xbp1KO")
names(listDEgenes)<-c("cDC1.Xbp1KO-cDC1.WT", "cDC1.Ire1KO-cDC1.WT","cDC1.Ire1KO-cDC1.Xbp1KO")

names(listDEgenes_moreStrict)<-names(listDEgenes)
names(listDEgenes_lessStrict)<-names(listDEgenes)

##### Get numbers of DE genes
lapply(listDEgenes,dim)
lapply(listDEgenes_lessStrict,dim)
lapply(listDEgenes_moreStrict,dim)
lapply(listallgenes,dim)


### Add geneSymbol in column (for the export)
listDEgenes<-lapply(listDEgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})
listDEgenes_lessStrict<-lapply(listDEgenes_lessStrict,function(x){dplyr::mutate(x,'gene'=rownames(x))})
listallgenes<-lapply(listallgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})

########################################
##### Write results
########################################

write.xlsx(listDEgenes, file = "results_Clint_2023/summary_DEgenes_clint_onlycDC1.xlsx")
write.xlsx(listDEgenes_lessStrict, file = "results_Clint_2023/summary_DEgenes_lessStrict_clint_onlycDC1.xlsx")
write.xlsx(listallgenes, file = "results_Clint_2023/summary_allgenes_clint_onlycDC1_2023.xls")

### Save results
saveRDS(listDEgenes, file="results_Clint_2023/Robjects/listDEgenes_clint_onlycDC1.rds")
saveRDS(listDEgenes_lessStrict, file="results_Clint_2023/Robjects/listDEgenes_lessStrict_clint_onlycDC1.rds")
saveRDS(listallgenes, file="results_Clint_2023/Robjects/listallgenes_clint_onlycDC1_2023.rds")

## Update August 2023: upload to GEO (without cDC2 samples!!!)
dir.create("results_Clint_2023/GEO")

# Check names
all(colnames(countData) == colnames(expTable))

write.table(countData, "results_Clint_2023/GEO/Bulk_miRNA_seq_raw_counts_cDC1s_Ire1_Xbp1_experiment.txt", sep="\t")
write.table(expTable, "results_Clint_2023/GEO/Bulk_miRNA_seq_normalized_counts_cDC1s_Ire1_Xbp1_experiment.txt", sep="\t")

################################################################################
########## GET ALL GENES DOWNSTREAM OF THE MIRNA's
################################################################################

########################################
##### Load data
########################################

# https://bioconductor.org/packages/release/bioc/vignettes/miRBaseConverter/inst/doc/miRBaseConverter-vignette.html
## Extra: check version MirBase
data(miRNATest)
miRNANames = miRNATest$miRNA_Name
version=checkMiRNAVersion(miRNANames, verbose = TRUE) ## v18

##### DE genes #####
DEmirs<-unique(c(DEgenesGroup1$miRNAname, DEgenesGroup2$miRNAname, DEgenesGroup3$miRNAname))

# ## Extra: Little test
# result1 = miRNA_NameToAccession(DEmirs,version = "v18")
# Accessions=result1$Accession
# Family_Info2=checkMiRNAFamily(Accessions)
# head(Family_Info2)
# result_sequence = getMiRNASequence(Accessions,targetVersion = "v18")
# head(result_sequence)

##### Info about miR family #####
family_info = read.table("documentation/miR_Family_Info.txt", sep="\t", stringsAsFactors = F, header = T)
family_info_mm = family_info %>% filter(`Species.ID` %in% c(10090))
family_info_mm[grep('mmu-miR-92',family_info_mm$MiRBase.ID),]

##### Genes downstream of (conserved) miR family #####
mirna = read.table("documentation/Conserved_Family_Info.txt", sep="\t", stringsAsFactors = F, header = T) %>% filter(`Species.ID` %in% c(10090))
mirna %>% group_by(`miR.Family`) %>% summarize(n=n()) %>% arrange(-n)

mirna$symbol = mirna$`Gene.Symbol`
mirna = mirna %>% drop_na()

##create list with all miRNA families and their targets
mirna = mirna %>% dplyr::rename(mir=`miR.Family`)
gsets = mirna %>% plyr::dlply(plyr::.(mir), function(x) x$symbol)

##### Genes downstream of (non-conserved) miR family #####
mirna_nc = read.table("documentation/Nonconserved_Family_Info.txt", sep="\t", stringsAsFactors = F, header = T) %>% filter(`Species.ID` %in% c(10090))
mirna_nc %>% group_by(`miR.Family`) %>% summarize(n=n()) %>% arrange(-n)

mirna_nc$symbol = mirna_nc$`Gene.Symbol`
mirna_nc = mirna_nc %>% drop_na()

##create list with all miRNA families and their targets
mirna_nc = mirna_nc %>% dplyr::rename(mir=`miR.Family`)
gsets_nc = mirna_nc %>% plyr::dlply(plyr::.(mir), function(x) x$symbol)

########################################
##### Get miR family
########################################

family_info_mm[grep('mmu-miR-3963',family_info_mm$MiRBase.ID),]
# miR-3963 (mmu-miR-3963)

head(gsets_nc$`miR-3963`)

family_info_mm[grep('mmu-miR-18a',family_info_mm$MiRBase.ID),]
# miR-18-5p (mmu-miR-18a-5p)
# miR-18-3p/7069-3p (mmu-miR-18a-3p)

head(gsets$`miR-18-5p`)
head(gsets_nc$`miR-18-3p/7069-3p`)

family_info_mm[grep('mmu-miR-10',family_info_mm$MiRBase.ID),]
# miR-10a-3p (mmu-miR-10a-3p)
# miR-10-5p (mmu-miR-10a-5p + mmu-miR-10b-5p)
# miR-10b-3p (mmu-miR-10b-3p)

head(gsets_nc$`miR-10a-3p`)
head(gsets$`miR-10-5p`)
head(gsets_nc$`miR-10b-3p`)

# family_info_mm[grep('mmu-miR-148a',family_info_mm$MiRBase.ID),]
# # miR-148-5p (mmu-miR-148a-5p)
# # miR-148-3p/152-3p (mmu-miR-148a-3p)
# 
# head(gsets_nc$`miR-148-5p`)
# head(gsets$`miR-148-3p/152-3p`)

family_info_mm[grep('mmu-miR-183',family_info_mm$MiRBase.ID),]
# miR-183-5p (mmu-miR-183-5p)
# miR-183-3p (mmu-miR-183-3p)
# miR-183-5p.2 (mmu-miR-183-5p.2)

head(gsets$`miR-183-5p`)
head(gsets_nc$`miR-183-3p`)
head(gsets$`miR-183-5p.2`)

family_info_mm[grep('mmu-miR-155',family_info_mm$MiRBase.ID),]
# miR-155-5p (mmu-miR-155-5p)
# miR-155-3p (mmu-miR-155-3p)

head(gsets$`miR-155-5p`)
head(gsets_nc$`miR-155-3p`)

family_info_mm[grep('mmu-miR-451',family_info_mm$MiRBase.ID),]
# miR-451a (mmu-miR-451a)
# miR-451b (mmu-miR-451b)

head(gsets$`miR-451a`)
head(gsets$`miR-451b`)

family_info_mm[grep('mmu-miR-455',family_info_mm$MiRBase.ID),]
# miR-455-5p (mmu-miR-455-5p)
# miR-455-3p.1 (mmu-miR-455-3p.1)
# miR-455-3p.2 (mmu-miR-455-3p.2)

head(gsets$`miR-455-5p`)
head(gsets$`miR-455-3p.1`)
head(gsets$`miR-455-3p.2`)

#######

family_info_mm[grep('mmu-miR-6935',family_info_mm$MiRBase.ID),]
# miR-705/6922-5p/6935-5p/6993-5p/7008-5p (mmu-miR-6935-5p)
# miR-6935-3p/7092-3p (mmu-miR-6935-3p)

head(gsets$`miR-6935-5p`)
head(gsets_nc$`miR-6935-3p`)

family_info_mm[grep('mmu-miR-703',family_info_mm$MiRBase.ID),]
# miR-703 (mmu-miR-703)

head(gsets$`miR-703`)

family_info_mm[grep('mmu-miR-18a',family_info_mm$MiRBase.ID),]
# miR-18-5p (mmu-miR-18a-5p)
# miR-18-3p/7069-3p (mmu-miR-18a-3p)

head(gsets$`miR-18-5p`)
head(gsets_nc$`miR-18-3p/7069-3p`)

family_info_mm[grep('mmu-miR-92',family_info_mm$MiRBase.ID),]
# miR-25-3p/32-5p/92-3p/363-3p/367-3p (mmu-miR-92a-3p + mmu-miR-92b-3p)
# miR-92b-5p (mmu-miR-92b-5p)
# miR-92a-2-5p (mmu-miR-92a-2-5p)
# miR-92a-1-5p (mmu-miR-92a-1-5p)

head(gsets$`miR-25-3p/32-5p/92-3p/363-3p/367-3p`)
head(gsets_nc$`miR-92b-5p`)
head(gsets_nc$`miR-92a-2-5p`)
head(gsets_nc$`miR-92a-1-5p`)


family_info_mm[grep('mmu-miR-10',family_info_mm$MiRBase.ID),]
# miR-10a-3p (mmu-miR-10a-3p)
# miR-10-5p (mmu-miR-10a-5p + mmu-miR-10b-5p)
# miR-10b-3p (mmu-miR-10b-3p)

head(gsets_nc$`miR-10a-3p`)
head(gsets$`miR-10-5p`)
head(gsets_nc$`miR-10b-3p`)


family_info_mm[grep('mmu-miR-155',family_info_mm$MiRBase.ID),]
# miR-155-5p (mmu-miR-155-5p)
# miR-155-3p (mmu-miR-155-3p)

head(gsets$`miR-155-5p`)
head(gsets_nc$`miR-155-3p`)

####################

myListGenes<-list("miR-3963" = gsets_nc$`miR-3963`,
                  "miR-18-3p/7069-3p" = gsets_nc$`miR-18-3p/7069-3p`,
                  "miR-10-5p" = gsets$`miR-10-5p`,
                  "miR-183-5p" = gsets$`miR-183-5p`,
                  "miR-155-5p" = gsets$`miR-155-5p`,
                  "miR-451a" = gsets$`miR-451a`,
                  "miR-455-3p.1" = gsets$`miR-455-3p.1`,
                  "miR-455-3p.2" = gsets$`miR-455-3p.2`,
                  # Second comparison
                  'miR-6935-3p' =gsets_nc$`miR-6935-3p`,
                  'miR-703' =gsets_nc$`miR-703`,
                  'miR-18a-3p' = gsets_nc$`miR-18a-3p`,
                  "miR-92a-1-5p" = gsets_nc$`miR-92a-1-5p`,
                  "miR-10-5p" = gsets$`miR-10-5p`,
                  "miR-155-5p" = gsets$`miR-155-5p`)

saveRDS(myListGenes, file="results_Clint_2023/Robjects/myListGenes.rds")

################################################################################
########## CHECK SOME KNOWN MIRNA's
################################################################################

########################################
##### Previous analysis
########################################

##### Starting from the 15 miRNAs #####
gsetidsoi<-readRDS(file="/home/clintdn/VIB/DATA/Sophie/RNA-seq_Simon/results/memeResults/gsetidsoi.rds")
gsetidsoi

family_info = read.table("documentation/miR_Family_Info.txt", sep="\t", stringsAsFactors = F, header = T)
family_info_oi = family_info %>% filter(`Species.ID` %in% c(10090)) %>% filter(`miR.family` %in% gsetidsoi)

allDEgenes<-unique(c(rownames(DEgenesGroup1), rownames(DEgenesGroup2), rownames(DEgenesGroup3)))
allDEgenes_lowS<-unique(c(rownames(DEgenesGroup1_lowS), rownames(DEgenesGroup2_lowS), rownames(DEgenesGroup3_lowS)))
allDEgenes_lessStrict<-unique(c(rownames(DEgenesGroup1_clint), rownames(DEgenesGroup2_clint), rownames(DEgenesGroup3_clint)))

intersect(allDEgenes, family_info_oi$MiRBase.Accession)
intersect(allDEgenes_lowS, family_info_oi$MiRBase.Accession)
intersect(allDEgenes_lessStrict, family_info_oi$MiRBase.Accession) #"MIMAT0000539" = 92a-3p


####################################################################################################
############################## ----------- TRIWISE PLOTS ------------ ##############################
####################################################################################################

##### All DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1_clint),rownames(DEgenesGroup2_clint),rownames(DEgenesGroup3_clint)))
length(allDEgenes)
##58

######################################################################
########## PREPARATION
######################################################################
wantedColors7<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red", 
                 nodiffSet2="limegreen",diffSet2="darkgreen", nodiffSet3="lightblue", diffSet3="darkblue", 
                 nodiffSet4="darkorchid1", diffSet4="darkmagenta", nodiffSet5="lightsalmon", diffSet5="darkorange", 
                 nodiffSet6="lightskyblue", diffSet6="turquoise4", nodiffSet7="bisque", diffSet7="burlywood4")
wantedColors<-c(nodiffall="gray55",diffall="indianred1")

##### Triwise plots
colsWT<-grep("Xbp1_WT|Ire1Xbp1_WT",colnames(expTable))
colsIre1KO<-grep("Ire1Xbp1_KO",colnames(expTable))
colsXpb1KO<-grep("^Xbp1_KO",colnames(expTable))

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

saveRDS(barycoords, "results_Clint_2023/Robjects/barycoords.rds")
saveRDS(Gdiffexp, "results_Clint_2023/Robjects/Gdiffexp.rds")

#######

## Check position DE miRNAs in WTvsDKO
## Choose list!!
barycoords_clint<-barycoords[rownames(DEgenesGroup2),] #new paper

barycoords_clint$miRNA<-rownames(barycoords_clint)

write.xlsx(barycoords_clint, "results_Clint_2023/triwiseResults/barycoords_clint_new_paper.xlsx") #new paper

p<-plotDotplot_byLM(theBarycoords, Gdiffexp=rownames(DEgenesGroup2), colorvalues=wantedColors, showlabels = F) #new paper

pdf(file="results_Clint_2023/triwiseResults/triwisePlot_DKO_vs_WT_new_paper.pdf",width =10, height = 8)
print(p)
dev.off()

## Truncate miRNA names
miRNAs_DKO_vs_WT<-gsub("mmu-","",DEgenesGroup2$miRNAname)
miRNAs_DKO_vs_WT[5]<-"miR-10-5p"
Genes_DKO_vs_WT<-myListGenes[miRNAs_DKO_vs_WT]

saveRDS(myListGenes[miRNAs_DKO_vs_WT], file="results_Clint_2023/Robjects/myListGenes_DKO_vs_WT.rds")

###########################################################################################################################################
###########################################################################################################################################

###Triwise plot
p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors, showlabels = F)
print(p)

p<-plotDotplot_byLM(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors, showlabels = F)
print(p)
ggsave(p,file="results_Clint_2023/triwiseResults/triwisePlot_new.png",dpi=300)


###Triwise plot advanced
overlap1and2<-intersect(rownames(DEgenesGroup1_lowS),rownames(DEgenesGroup2_lowS))
overlap1and3<-intersect(rownames(DEgenesGroup1_lowS),rownames(DEgenesGroup3_lowS))
overlap2and3<-intersect(rownames(DEgenesGroup2_lowS),rownames(DEgenesGroup3_lowS))

only1<-setdiff(rownames(DEgenesGroup1_lowS),c(overlap1and2,overlap1and3))
only2<-setdiff(rownames(DEgenesGroup2_lowS),c(overlap1and2,overlap2and3))
only3<-setdiff(rownames(DEgenesGroup3_lowS),c(overlap1and3,overlap2and3))

p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F, Goi=list(Set1=only1,Set2=only2,Set3=only3,Set4=overlap1and2,
                                                                            Set5=overlap1and3,Set6=overlap2and3),
               colorvalues=wantedColors7) +
  theme(legend.position = "right")
print(p)
ggsave(p,file="results_Clint_2023/triwiseResults/triwisePlot_advanced.png",dpi=300)


###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)
ggsave(p,file="results_Clint_2023/triwiseResults/rosePlot_new.png",dpi=300)


p<-plotRoseplot(theBarycoords, Gdiffexp=rownames(DEgenesGroup2_clint), showlabels = F)
print(p)
ggsave(p,file="results_Clint_2023/triwiseResults/rosePlot_new_DE_genes_DKO_vs_WT.png",dpi=300)

## Extra Victor (22/03/23)
rownames(expTable_mean)
expTable_mean_Victor<-expTable_mean
rownames(expTable_mean_Victor) = miRNA_AccessionToName(rownames(expTable_mean_Victor))[,2]

mIRS_Victor<-c("mmu-miR-200b-3p", "mmu-miR-200a-3p", "mmu-miR-200c-3p", "mmu-miR-200a-5p",
               "mmu-miR-34a-5p", "mmu-miR-34a-3p")
###Interactive plot
p<-interactiveDotplot(expTable_mean_Victor, Gdiffexp=rownames(DEgenesGroup2_clint), Goi = mIRS_Victor, plotLocalEnrichment=FALSE)
print(p) 
saveWidget(p,file="/home/clintdn/VIB/DATA/Sophie/miRNASeq/results_Clint_2023/triwiseResults/interactiveTriwisePlot_Victor_22032023.html") ##needs full path!

## Extra Victor thesis (02/08/23)
rownames(expTable_mean)
expTable_mean_Victor<-expTable_mean
rownames(expTable_mean_Victor) = miRNA_AccessionToName(rownames(expTable_mean_Victor))[,2]
mIRS_Victor<-DEgenesGroup2$miRNAname

###Interactive plot
p<-interactiveDotplot(expTable_mean_Victor, Gdiffexp=rownames(DEgenesGroup2), Goi = mIRS_Victor, plotLocalEnrichment=FALSE)
print(p) 
saveWidget(p,file="/home/clintdn/VIB/DATA/Sophie/miRNASeq/results_Clint_2023/triwiseResults/interactiveTriwisePlot_Victor_02082023.html") ##needs full path!
