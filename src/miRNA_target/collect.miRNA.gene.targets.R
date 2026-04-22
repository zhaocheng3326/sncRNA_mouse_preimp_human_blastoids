#' ---
#' title: collect miRNA targets
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
  #library(ff)
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
load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)


targetscan.pred <- read.delim("/home/chenzh/Genome_new/targetscan/old/mmu.miRNA.target.long.txt",head=F,stringsAsFactors = F,sep = "\t") %>% tbl_df() %>% select(V2,V3) %>% setNames(c("gene","miRNA"))  %>% group_by(miRNA,gene) %>% summarise(nTargets=n())%>% ungroup()%>% mutate(gene=paste(toupper(substring(gene, 1, 1)), tolower(substring(gene, 2)), sep=""))
miRDB.pred <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/miRDB/miRDB.mmu.miRNA.target.gene.out",,head=T,stringsAsFactors = F,sep = "\t") %>% tbl_df() %>% mutate(gene=Gene) %>% select(gene,miRNA)%>% group_by(miRNA,gene) %>% summarise(nTargets=n()) %>% ungroup() 

target.pairs <- miRDB.pred %>% anti_join(targetscan.pred,by = c("miRNA", "gene")) %>% bind_rows(anti_join(targetscan.pred, miRDB.pred,by = c("miRNA", "gene"))) %>% bind_rows(miRDB.pred %>% inner_join(targetscan.pred,by = c("miRNA", "gene")) %>% mutate(nTargets=ifelse(nTargets.x > nTargets.y,nTargets.x,nTargets.y)) %>% select(-c(nTargets.x,nTargets.y))) %>% group_by(miRNA,gene) %>% top_n(1,nTargets)

target.pairs <- target.pairs %>% filter(miRNA %in% miRNA.bed$V7)%>%  ungroup() %>% filter(gene %in% Gene.discrp$V6)
saveRDS(target.pairs,paste0("tmp_data/",TD,"/miRNA.gene.target.rds"))
