#' ---
#' title: cor for gene and miRNA for coseq
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
  library(scran)
  #library(WGCNA)
  #library(hdWGCNA)
  #library(UCell)
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

TF_miRNA <- readRDS(paste0("tmp_data/",TD,"/TF_miRNA.rds")) ## TF regulating miRNA
target.pairs <- readRDS(paste0("tmp_data/",TD,"/miRNA.gene.target.rds")) #miRNA targetting gene
# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=recode(EML,"prelineage"="L8CM")) %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% select(cell,EML) %>% filter(EML!="unknown") %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the full cell annotation
meta.filter <- meta.filter  %>% inner_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap %>% filter(EML!="unknown"),by="cell")  %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))

sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% pull(cell) 
blast_sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% filter(EML %in% c("ICM","TE")) %>% pull(cell) 
#' miRNA exp
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))[,sel_cells]
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))


if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.rds"))) {
  miRNA.gene.cor.out <- readRDS(file=paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.rds"))
  miRNA.gene.blast.cor.out <- readRDS(file=paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.rds"))
}else{
  #' calculate correlation between miRNA and gene modules
  left_data <- as.data.frame(coseq.miRNA.exp[,sel_cells])
  right_data <- as.data.frame(coseq.gene.exp[,sel_cells])
  left_right <- target.pairs%>% select(miRNA,gene) %>% setNames(c("left","right"))
  miRNA_gene.cor <- FunPairedCor(left_data ,right_data,left_right,method="spearman") %>% ungroup() %>% rename(miRNA=left,gene=right)
  
  
  left_data <- as.data.frame(coseq.miRNA.exp[,blast_sel_cells])
  right_data <- as.data.frame(coseq.gene.exp[,blast_sel_cells])
  left_right <- target.pairs%>% select(miRNA,gene) %>% setNames(c("left","right"))
  miRNA_gene.blast.cor <- FunPairedCor(left_data ,right_data,left_right,method="spearman") %>% ungroup() %>% rename(miRNA=left,gene=right)
  
  
  
  saveRDS(miRNA_gene.cor,file=paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.rds"))
  saveRDS(miRNA_gene.blast.cor,file=paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.rds"))
}


