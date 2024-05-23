# Spleen_cDC1_Homeostatic_Maturation_in_Ire1_Xbp1_DKO

## About

Conventional type I dendritic cells (cDC1s) are characterized by their unique capacity to engulf and cross-present antigens from dead cells, in contrast to the closely related cDC2s. Recently, the engulfment of dead cells by cDC1s was found to initiate cDC1 homeostatic maturation and induction of cholesterol efflux pathways by inducing LXRb activity. cDC1s are uniquely characterized by high basal activity of IRE1, a sensor of the unfolded protein response (UPR). IRE1 activity in cDC1s is unconventional because it does not induce a typical UPR gene signature and hence the role of IRE1 in cDC1s remains enigmatic. Here we show that IRE1 activity is critical for the survival of homeostatically matured cDC1s, not cDC2s, after engulfment of apoptotic cells. Mechanistically, we found that IRE1 becomes activated by the influx of dead-cell derived cholesterol, explaining its cDC1 subset specific activation. Following activation, IRE1’s endonuclease activity degrades miRNA-92a in a process termed Regulated IRE1 dependent decay (RIDD). miRNA-92a targets Abcg1 mRNA which is an essential cholesterol efflux transporter in cDC1s. Thus IRE1, together with LXRß which drives the expression of Abcg1 mRNA, coordinates cholesterol efflux in cDC1s after efferocytosis. Loss of IRE1 leads to cholesterol toxicity and cell death of mature cDC1s after engulfment of dead cells. This could be rescued by blocking miRNA synthesis or by enforcing cholesterol efflux with reconstituted high density lipoprotein treatment. In conclusion, these data highlight the central role of the endoplasmic reticulum and IRE1 as a sensor of cholesterol influx, extending IRE1’s function beyond its canonical role in protein folding. Furthermore, they underscore the tight control of lipid metabolism during cDC1 maturation, uncovering a second pathway to coordinate cholesterol efflux.


## Overview scripts

Here's an overview of the various R scripts used in processing the Bulk RNA-seq, miRNA-seq and CITE-Seq data in the manuscript Bosteels et al.:
- 1.script_BulkRNAseq_cDC1s_triwise.R: Limma-edgeR workflow for running differential expression analysis and triwise analysis on the cDC1 Bulk RNA-seq data
- 2.script_BulkRNAseq_cDC1s_findMiRNAfamilies.R: Script for analyzing the miRNA target genes and subsequent motif analysis in the cDC1 Bulk RNA-seq data 
- 3.script_BulkRNAseq_cDC2s.R: Limma-edgeR workflow for running differential expression analysis on the cDC2 Bulk RNA-seq data
- 4.script_CITEseq_BareBones_RNA-ADT_HPCscript.R: Standard pipeline used on the High Performance Cluster for analyzing the RNA and ADT assays of the CITE-seq objects
- 5.script_CITEseq_merge_part1.R: Script for basic analysis of the cDC1 subset of merged WT and DKO CITE-seq data 
- 6.script_CITEseq_merge_part2.R: Script for downstream analysis of the cDC1 subset of merged WT and DKO CITE-seq data
- 7.script_CITEseq_merge_muscat.R: Script for muscat Differential State analysis between DKO and WT in the merged CITE-seq data
- 8.script_miRNAseq.R: Limma-edgeR workflow for running differential expression analysis on the cDC1 miRNA-seq data

## Citation

Victor Bosteels et al., The UPR sensor IRE1a is essential for homeostatic dendritic cell maturation.
