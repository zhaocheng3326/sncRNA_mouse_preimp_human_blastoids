#' ---
#' title: "QC". 
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
#tRNA.nExpG.cutoff <- 100

sel.type <- c("miRNA","snoRNA","snRNA","tRNA","rRNA","piRNA") #,"os_piRNA"

# ##Loading transcript annotation
if (file.exists(paste0("tmp_data/",TD,"/Gene.meta.Rdata"))) {
  load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)
}

td="merge_batch_soft_link"
### Loading expression matrix and mapping information
if (file.exists(paste0("tmp_data/",td,"/all.raw.data.Rdata"))) {
  load(paste0("tmp_data/",td,"/all.raw.data.Rdata"),verbose=T)
}


#' QC 
#meta.filter <- meta.all %>% filter(devTime!="sperm") %>% filter(Total_reads >Total_reads.cutoff & mt.perc < mt.perc.cutoff  & miRNA.nExpG >miRNA.nExpG.cutoff)  %>% bind_rows(meta.all %>% filter(devTime=="sperm") %>% filter(Total_reads >Total_reads.cutoff  & miRNA.nExpG >miRNA.nExpG.cutoff)) %>% filter(!cell %in% em_H9.minor.cells) # & tRNA.nExpG > tRNA.nExpG.cutoff)
meta.filter <- meta.all %>% filter(devTime!="sperm") %>% filter(Total_reads >Total_reads.cutoff & mt.perc < mt.perc.cutoff  & miRNA.nExpG >miRNA.nExpG.cutoff)  %>% bind_rows(meta.all %>% filter(devTime=="sperm") %>% filter(Total_reads >Total_reads.cutoff  & miRNA.nExpG >miRNA.nExpG.cutoff)) %>% filter(!embryo %in% c("em_H9")) # & tRNA.nExpG > tRNA.nExpG.cutoff)

#' remove the low quality embryo and cells
meta.filter <- meta.filter %>%filter(!embryo %in% batch5.lq.embyro.ID)%>%filter(!embryo %in% split2.lq.embyro.ID )%>%filter(!cell %in% lq.small.cells )

#' only keep the 2C and 4C cells for Split2
meta.filter <- meta.filter %>% filter(batch!="Split2") %>% bind_rows(meta.filter %>% filter(devTime %in% c("2C","4C") & batch=="Split2"))



counts.filter <- counts
counts.filter$mature <- counts$mature %>% select(ID,meta.filter$cell)

meta.filter$chrX.perc <- colSums(counts.filter$mature[chrX.miRNA,meta.filter$cell])/meta.filter$miRNA.UMI
meta.filter$chrX.unique.perc <- colSums(counts.filter$mature[chrX.miRNA.unique,meta.filter$cell])/meta.filter$miRNA.UMI
meta.filter$chrY.perc <- NA
meta.filter$chrY.unique.perc <- NA



rewrite=FALSE
if (rewrite) {
  saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/small.meta.filter.rds"))
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/small.counts.filter.rds"))
}
  


