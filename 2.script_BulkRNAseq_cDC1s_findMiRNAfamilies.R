## Script for analyzing Xbp1/Ire1KO cDC1 bulk RNA-seq project data
## Check target genes of miRNA-families on triwise and check common motifs

library("tidyverse")
library("triwise")
library("magrittr")
library("mirbase.db")
library('Biostrings')
library('PWMEnrich')
library('xml2')

################################################################################
######### LOAD DATA
################################################################################

getwd()

#####  Triwise data
barycoords <- readRDS("results/barycoords.rds")
Gdiffexp <- readRDS("results/Gdiffexp.rds")


#####  miRNA data
mirna = read_tsv("rawData/Conserved_Family_Conserved_Targets_Info.txt.zip")# %>% filter(`Species ID`==10090)
mirna %>% group_by(`miR Family`) %>% summarize(n=n()) %>% arrange(-n)
# `miR Family`     n
# <chr> <int>
# 1    miR-30-5p  9640
# 2     miR-9-5p  8993
# 3    miR-27-3p  8950
# 4   miR-181-5p  8155
# 5    miR-29-3p  7437
# 6    miR-19-3p  7284
# 7      miR-182  7244
# 8    miR-23-3p  7019
# 9   miR-340-5p  6584
# 10   miR-218-5p  5889

################################################################################
######### PROCESS DATA
################################################################################

### Convert Ensembl into Entrez
entrez2ensembl = as.list(org.Mm.eg.db::org.Mm.egENSEMBL)
entrez2ensembl = entrez2ensembl %>% keep(~!is.na(.[[1]])) %>% map_chr(~.[[1]])
ensembl2entrez = setNames(names(entrez2ensembl), entrez2ensembl)

mirna$ensembl = gsub("(.*)\\..*", "\\1", mirna$`Gene ID`)
mirna$entrez = ensembl2entrez[mirna$ensembl]
mirna$symbol = mirna$`Gene Symbol`
mirna = mirna %>% drop_na()

### Create list with all miRNA families and their targets
mirna = mirna %>% dplyr::rename(mir=`miR Family`)
gsets = mirna %>% plyr::dlply(plyr::.(mir), function(x) x$symbol)
gsets = sapply(gsets, function(gset) intersect(rownames(barycoords), unique(as.character(gset))))

### Histogram of number of targets per miRNA family
map_int(gsets, length) %>% hist


################################################################################
######### PLOT miRNA FAMILIES ON TRIWISE PLOT
################################################################################

plotDotplot(barycoords, Gdiffexp=Gdiffexp)
plotRoseplot(barycoords, Gdiffexp = Gdiffexp)


### Take all miRNA targets; check if there are enriched in a certain direction of the triwise plot
scores = testUnidirectionality(barycoords, gsets, Gdiffexp=Gdiffexp, statistic="angle")
scores<-scores[order(scores$qval),]
head(scores,10)

cutOff_qval<-0.15
sigScores<-scores %>% filter(qval <= cutOff_qval)
sigScores

### Plot scores
plotPvalplot(scores, attributes(barycoords)$conditions)
plotPvalplot(sigScores, attributes(barycoords)$conditions)

### Plot miRNA families
plots = lapply(scores %>% head(8) %>% .$gsetid, function(gsetid) {
  plotDotplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::ggtitle(gsetid %>% strwrap(40) %>% paste(collapse="\n")) + 
    ggplot2::theme(axis.text=ggplot2::element_text(size=14), plot.title=ggplot2::element_text(size=8,face="bold"))
})
cowplot::plot_grid(plotlist=plots, ncol=4)


plots = lapply(scores %>% head(8) %>% .$gsetid, function(gsetid) {
  plotRoseplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::ggtitle(gsetid %>% strwrap(40) %>% paste(collapse="\n")) + 
    ggplot2::theme(axis.text=ggplot2::element_text(size=14), plot.title=ggplot2::element_text(size=8,face="bold"))
})
cowplot::plot_grid(plotlist=plots, ncol=4)


# ################################################################################
# ######### WORK ON CERTAIN CORNER OF TRIWISE PLOT
# ################################################################################
# 
# ### Select bottom corner
# tmp<-barycoords[barycoords$angle < -0.3,]
# barycoordsSlice<-tmp[tmp$angle > -2.5,]
# p<-plotDotplot(barycoords, Gdiffexp=Gdiffexp, showlabels = F, Goi=list(Set1=rownames(barycoordsSlice)))
# print(p)
# 
# ### Select right top corner
# tmp<-barycoords[barycoords$angle < 1.7,]
# barycoordsSlice<-tmp[tmp$angle > 0.3,]
# p<-plotDotplot(barycoords, Gdiffexp=Gdiffexp, showlabels = F, Goi=list(Set1=rownames(barycoordsSlice)))
# print(p)
# 
# 
# 
# scores <- testEnrichment(intersect(rownames(barycoordsSlice), Gdiffexp), gsets %>% map(~intersect(., Gdiffexp)), Gdiffexp)
# scores %>% arrange(qval) %>% head(.,10)
# 
# cutOff_qval<-0.01
# sigScores<-scores %>% filter(qval <= cutOff_qval)
# sigScores

################################################################################
######### GET miRNA families OF INTEREST
################################################################################

### Get the miRNA families of interest
gsetidsoi = as.character(scores %>% filter(qval < cutOff_qval) %>% .$gsetid)
scores[scores$gsetid %in% gsetidsoi,] %>% arrange(qval)
saveRDS(gsetidsoi,"results/memeResults/gsetidsoi.rds")

### Create triwise plot and rose plot of only the miRNA families of interest
plots = lapply(gsetidsoi, function(gsetid) {
  plotDotplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::ggtitle(gsetid %>% strwrap(40) %>% paste(collapse="\n")) + 
    ggplot2::theme(axis.text=ggplot2::element_text(size=14), plot.title=ggplot2::element_text(size=8,face="bold"))
})
cowplot::plot_grid(plotlist=plots, ncol=4)

plots = lapply(gsetidsoi, function(gsetid) {
  plotRoseplot(barycoords, Gdiffexp, Goi=gsets[[gsetid]], showlabels = F) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::ggtitle(gsetid %>% strwrap(40) %>% paste(collapse="\n")) + 
    ggplot2::theme(axis.text=ggplot2::element_text(size=14), plot.title=ggplot2::element_text(size=8,face="bold"))
})
cowplot::plot_grid(plotlist=plots, ncol=4)

### Create heatmap to show overlap in genes behind the miRNA families of interest
jaccard = function(x1, x2) sapply(1:length(x1), function(i) length(intersect(gsets[[x1[[i]]]], gsets[[x2[[i]]]]))/length(union(gsets[[x1[[i]]]], gsets[[x2[[i]]]])))
symmetric = function(x) {x[lower.tri(x)] = t(x)[lower.tri(x)]; x}
combinations = combn(gsetidsoi, 2) %>% t %>% as.data.frame %>% mutate(jaccard=jaccard(V1, V2))
combinations %>% reshape2::acast(V1 ~ V2, value.var = "jaccard") %>% symmetric() %>% pheatmap::pheatmap()

### Take all genes behind the miRNA families of interest and create rose plot of it
Goi <- unlist(gsets[gsetidsoi])
subbarycoords = barycoords[unique(Goi),]
subbarycoords$number = table(Goi)[rownames(subbarycoords)]
subbarycoords %<>% tibble::rownames_to_column("gene") %>% mutate(diffexp=gene %in% Gdiffexp)

plotRoseplot(subbarycoords %>% filter(diffexp), Coi=c("", "", ""))

### Only take the genes that are behind all miRNA families or all miRNA families execpt one and so one ...
Goitable = tibble(gene=Goi, mir=names(Goi))
subbarycoords$mirnas = names(Goi) %>% tapply(Goi, function(x) list(x)) %>% .[subbarycoords$gene]
subbarycoords %<>% mutate(n=map_int(mirnas, ~length(.)))
subbarycoords %<>% mutate(mirnas_text=map_chr(mirnas, ~paste0(., collapse=" ")))

test<-subbarycoords %>% dplyr::select(-mirnas) %>% arrange(-number) %>% filter(diffexp)
head(test)

maxNumberToCheck<-12

test<-map(1:maxNumberToCheck, function(numberoi) {
  plotRoseplot(subbarycoords %>% filter(diffexp) %>% filter(number>=numberoi), Coi=c("", "", ""))
}) %>% cowplot::plot_grid(plotlist=.)
print(test)
ggsave(test,file="results/memeResults/triwisePlots_about_miRNAs/piePlot_15miRNA_numberOfGenes.png", dpi=300)

#Each plot separatly
for(i in 1:maxNumberToCheck){
  test<-plotRoseplot(subbarycoords %>% filter(diffexp) %>% filter(number>=i), Coi=c("", "", ""))
  fileName<-paste0("results/memeResults/triwisePlots_about_miRNAs/piePlot_15miRNA_numberOfGenes_N",i,".png")
  ggsave(test,file=fileName, dpi=300)
}

#Write genes to Excel file
xlsx::write.xlsx(test, "results/memeResults/list_15miRNA_withGenes.xls", "targetGenes")

################################################################################
######### CHECK MOTIF IN miRNAs OF INTEREST
################################################################################
family_info = read_tsv("rawData/miR_Family_Info.txt")# %>% filter(`Species ID`==10090)

hairpins = mirbaseHAIRPIN %>% as.list()
sequences = mirbaseSEQUENCE %>% as.list()

family_info_oi = family_info %>% filter(`Species ID` %in% c(9606, 10090)) %>% filter(`miR family` %in% gsetidsoi)

mirbaseids = family_info_oi$`MiRBase ID` %>% tolower()
mirbaseids_alternative = mirbaseids %>% gsub("(.*)-.*$", "\\1", .) %>% tolower() # cut off 3p/5p stuff
mirbaseids = ifelse(mirbaseids %in% names(hairpins), mirbaseids, mirbaseids_alternative) %>% unique %>% keep(~. %in% names(hairpins))
#hairpins[mirbaseids] %>% map(~cat(., "\n"))
sequences[mirbaseids] %>% unlist() %>% stringr::str_count("CUGCAG")

sequences %>% unlist() %>% stringr::str_count("CUGCAG") %>% hist


################################################################################
######### CHECK PWM MOTIF IN miRNAs OF INTEREST
################################################################################

### Get sequences for the background
rnasequences = Biostrings::RNAStringSet(sequences %>% as.character)
names(rnasequences) = names(sequences)
dnasequences = Biostrings::DNAStringSet(sequences %>% as.character %>% gsub("U", "T", .))
names(dnasequences) = names(sequences)

### Create function to easily print
bold_hairpin = function(hairpin, motif_position, motif_width) {
  hairpin_split = strsplit(hairpin, "\n")[[1]]
  newhairpin = list("", hairpin_split[[3]], "")
  for (i in 1:max(nchar(hairpin_split[[1]]), nchar(hairpin_split[[2]]))) {
    if((hairpin_split[[1]] %>% substr(i, i)) %in% c("A", "U", "G", "C", "a", "c", "u", "g", "-")) {
      newhairpin[[1]] = paste0(newhairpin[[1]], hairpin_split[[1]] %>% substr(i, i))
    } else {
      newhairpin[[1]] = paste0(newhairpin[[1]], hairpin_split[[2]] %>% substr(i, i))
    }
    if((hairpin_split[[4]] %>% substr(i, i)) %in% c("A", "U", "G", "C", "a", "c", "u", "g", "-")) {
      newhairpin[[3]] = paste0(newhairpin[[3]], hairpin_split[[4]] %>% substr(i, i))
    } else {
      newhairpin[[3]] = paste0(newhairpin[[3]], hairpin_split[[5]] %>% substr(i, i))
    }
  }
  newhairpin[[3]] = reverse(newhairpin[[3]])
  newhairpin = paste(newhairpin, collapse="\n")
  
  
  hairpin_split = strsplit(newhairpin, "")[[1]]
  
  curposition = 0
  for (i in 1:length(hairpin_split)) {
    if(hairpin_split[[i]] %in% c("A", "U", "G", "C", "a", "c", "u", "g")) {
      curposition = curposition + 1
      if(curposition > motif_position && curposition < motif_position + motif_width + 1) {
        hairpin_split[[i]] = paste0("<b>", hairpin_split[[i]], "</b>")
      } else {
        hairpin_split[[i]] = paste0("", hairpin_split[[i]], "")
      }
    }
  }
  hairpin_split = paste0(hairpin_split, collapse="") %>% strsplit(., "\n") %>% .[[1]]
  hairpin_split[[3]] = reverse(hairpin_split[[3]]) %>% gsub(">b/<", "<b>", .) %>% gsub(">b<", "</b>", .)
  hairpin_split %>% paste0(collapse="\n")# %>% cat()
}

### Create PWM motif
motifs_ire1a = c("CUGCAG","CUGCCG","CUGCAG","CUGCAG","CAGCAG","CUGCAG","CUGCAA","CUGCAG","CUGCAG","CUGCAG","CCGCAG","CUGCAG") %>% gsub("U", "T", .) %>% DNAStringSet()
names(motifs_ire1a) = 1:length(motifs_ire1a)
pfm_ire1a = consensusMatrix(motifs_ire1a)[1:4,]
pwm_ire1a = pfm_ire1a %>% PWM()
print(pwm_ire1a)
# [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
# A 0.000 0.070 0.000 0.000 0.165 0.070
# C 0.169 0.070 0.000 0.169 0.070 0.000
# G 0.000 0.000 0.169 0.000 0.000 0.165
# T 0.000 0.161 0.000 0.000 0.000 0.000

### Check if motif is present
bg = makePWMEmpiricalBackground(dnasequences, list(pfm_ire1a, pfm_ire1a))
motifenrichment = motifEnrichment(DNAStringSet(sequences[mirbaseids] %>% as.character %>% gsub("U", "T", .)), bg)

groupReport(motifenrichment, top=1)$p.value
# 0.535 0.535

##print result
counter<-0
for(i in 1:length(mirbaseids)) {
  mirbaseid = mirbaseids[[i]]
  starts = matchPWM(pwm_ire1a, dnasequences[[mirbaseid]], min.score="80%") %>% start
  hairpin = hairpins[[mirbaseid]]

  print(starts)
  if(length(starts) > 0) {
    counter<-counter+1
    cat(mirbaseid)
    cat("\n")
    cat("\n")
    bold_hairpin(hairpin, starts[[1]], 6) %>% cat
    cat("\n")
    cat("\n")
  }
}
paste0("Number of miRNA's that contain the motif: ",counter)

# "Number of miRNA's that contain the motif: 4"

################################################################################
######### DE NOVO MOTIF DETECTION (MEME)
################################################################################

writeXStringSet(rnasequences, "results/memeResults/bg.fasta")
writeXStringSet(rnasequences[mirbaseids], "results/memeResults/test.fasta")

# COMMANDS:
# export PATH=$HOME/meme/bin:$PATH 
# psp-gen -pos test.fasta -neg bg.fasta > bg.psp
# meme test.fasta -rna -oc . -nostatus -time 17958 -maxsize 60000 -mod zoops -nmotifs 3 -minw 4 -maxw 10 -psp bg.psp

### Read results
xml = xml2::read_xml("results/memeResults/meme.xml")
xml_motif = xml %>% xml_find_first(".//motif[@id='motif_1']") 
motif_width = xml_motif %>% xml_attr("width") %>% as.numeric
xml_sites = xml_motif %>% xml_find_all(".//contributing_site")
motif_sites = xml_sites %>% map(function(xml_site) {
  data.frame(
    mirna=mirbaseids[(gsub(".*_([0-9]*)$", "\\1", xml_attr(xml_site, "sequence_id")) %>% as.numeric()) + 1],
    position=as.numeric(xml_attr(xml_site, "position"))
  )}) %>% bind_rows()


### Print results
cat("<pre>\n")
for(i in 1:nrow(motif_sites)) {
  row = motif_sites[i,]
  hairpin = hairpins[[row$mirna]]
  motif_position = row$position
  cat(row$mirna)
  cat("\n")
  cat("\n")
  bold_hairpin(hairpin, motif_position, motif_width) %>% cat
  cat("\n")
  cat("\n")
}
cat("</pre>")



motif_sites$sequence = lapply(1:nrow(motif_sites), function(i) {
  row = motif_sites[i,]
  sequence = rnasequences[[row$mirna]]
  motif_position = row$position
  
  substr(sequence, motif_position, motif_position + motif_width) %>% as.character() %>% unlist
})
motif_sites

### Create logo
motifs <- motifs_ire1a %>% gsub("T", "U", .) %>% as.character
motifs

foundMotifs<-as.matrix(motif_sites$sequence)
foundMotifs

### Create logo via http://weblogo.threeplusone.com/
