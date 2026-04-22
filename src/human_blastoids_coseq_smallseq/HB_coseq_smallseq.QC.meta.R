#' ---
#' title: "QC for HB_coseq small"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
rm(list=ls())
rewrite=FALSE

#' check whether in local computer
if (grepl("KI-",Sys.info()['nodename'])) {
  print("local computer")
  source("/Users/cheng.zhao/chzhao_bioinfo/PC/SnkM/SgCell.R")
  base_dir <- "/Users/cheng.zhao/Documents"
} else {
  print("On server")
  condaENV <- "/home/chenzh/miniconda3/envs/R4.3" 
  print(condaENV)
  #LBpath <- paste0(condaENV ,"/lib/R/library")
  #.libPaths(LBpath)
  base_dir="/home/chenzh"
}



suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  #library(scran)
  #library(batchelor)
  #library(SeuratWrappers)
  #library(slingshot)
  #library(tradeSeq)
  #library(ComplexHeatmap)
})


set.seed(1)
DIR <- paste0(base_dir,"/My_project/mouse_smallseq_preimp")
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)

## define the EM brief structure data
#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
#source("src/figures.setting.R")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="July_2025"

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

mt.perc.cutoff <- 0.1
Total_reads.cutoff <- 0.35*10^6
miRNA.nExpG.cutoff <- 100
tRNA.nExpG.cutoff <- 100

sel.type <- c("miRNA","snoRNA","snRNA","tRNA","rRNA","piRNA") #,"os_piRNA"

#' #### loading human gene meta
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)



td="HB_smallseq"
### Loading expression matrix and mapping information

load(paste0("tmp_data/",td,"/HB_coseq_smallseq.all.raw.data.Rdata"),verbose=T)

#' QC 
meta.filter <- meta.all %>% filter(Total_reads >Total_reads.cutoff & mt.perc < mt.perc.cutoff  & miRNA.nExpG >miRNA.nExpG.cutoff)  

#' exclude the blastoids 4-8 in smallseq
meta.filter <- meta.filter %>% filter( !(batch=="HB.smallseq" & embryo %in% c("blastoids4","blastoids5","blastoids6","blastoids7","blastoids8")))

counts.filter <- counts
counts.filter$mature <- counts$mature %>% select(ID,meta.filter$cell)
#counts.filter$RL.dis <- counts$RL.dis %>% select(type,len,meta.filter$cell)
#counts.filter$isform <- counts$isform %>% select(IS5p,IS3p,type,meta.filter$cell)
#counts.filter$modication <- counts$modication %>% select(combp,faT3bp,type,meta.filter$cell)
#counts.filter$bin.win <- counts$bin.win %>% select(ID,meta.filter$cell)


meta.filter$chrX.perc <- colSums(counts.filter$mature[counts.filter$mature$ID %in% chrX.miRNA,meta.filter$cell])/meta.filter$miRNA.UMI
meta.filter$chrX.unique.perc <- colSums(counts.filter$mature[counts.filter$mature$ID %in% chrX.miRNA.unique,meta.filter$cell])/meta.filter$miRNA.UMI
meta.filter$chrY.perc <- colSums(counts.filter$mature[counts.filter$mature$ID %in% chrY.miRNA,meta.filter$cell])/meta.filter$miRNA.UMI
meta.filter$chrY.unique.perc <- colSums(counts.filter$mature[counts.filter$mature$ID %in% chrY.miRNA.unique,meta.filter$cell])/meta.filter$miRNA.UMI

rewrite=FALSE
if (rewrite) {
  saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.meta.filter.rds"))
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))
}
  


