#' ---
#' title: "QC"
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
  condaENV <- "/home/chenzh/miniconda3/envs/R4.0" 
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

## 
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


qc.nGene.min <- 1000
qc.mt.perc <- 0.25


load(paste0("tmp_data/","/mouse_smt2seq","/mouse_coseq_smt2.all.counts.meta.Rdata"),verbose = T)
load(paste0("tmp_data/","/mouse_smt2seq","/mouse_coseq_smt2.gene.meta.Rdata"),verbose=T)

meta.all <- meta.all %>% filter(batch=="Split1") %>% bind_rows(meta.all %>% filter(batch=="Split2" & devTime %in% c("2C","4C")))



meta.filter <- meta.all %>% filter(  mt.perc < qc.mt.perc & nGene > qc.nGene.min )
counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter$cell]
  
saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/mouse_coseq_smt2.counts.filter.rds"))
saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/mouse_coseq_smt2.meta.filter.rds"))


#'
#chrY.genes <- gtf.anno %>% filter(V1=="chrY") %>% pull(V6) %>% unique()
temp1 <- read.delim("tmp_data/mouse_smt2seq/merge.rsem_rpkm.csv",sep=",",head=T,row.names = 1)
temp2 <-  read.delim("tmp_data/mouse_smt2seq_zygote//merge.rsem_rpkm.csv",sep=",",head=T,row.names = 1)

table(rownames(temp1)==rownames(temp2))

RSEM.exp <- ((temp1 %>% cbind(temp2)) %>% as.data.frame()) [,meta.filter$cell]
chrY.ENS.genes <- rownames(gtf.anno %>% filter(V1=="chrY")) %>% intersect(rownames(RSEM.exp)) %>% setdiff( rownames(gtf.anno %>% filter(V1!="chrY")))
chrY.genes <- gtf.anno[chrY.ENS.genes,"V6"] %>% unique()
# temp.para.list <- list()
# for (a in c(25,30,35,40,45,50,60,70)) {
#   for (b in c(75,80,90,100,120,125)) {
#     temp.para.list[[paste(a,b)]] <- meta.filter.mod %>% mutate(sex=ifelse(chrY.sum < a, "F","Undef")) %>% mutate(sex=ifelse(chrY.sum > b, "M",sex))  %>% group_by(embryo,sex) %>% summarise(nC=n()) %>% group_by(embryo) %>% mutate(prop=nC/sum(nC)) %>% group_by(embryo) %>% top_n(1,prop) %>% mutate(prop=1-prop) %>% pull(prop) %>% mean()
#   }
# }
# which.min(temp.para.list)
#' method one
exp.min <- 25
exp.max <- 75
exp.prop <- 0.5
meta.filter.mod <- meta.filter
meta.filter.mod$chrY.RSEM.sum <- colSums(RSEM.exp[chrY.ENS.genes,meta.filter$cell])
meta.filter.mod$chrY.ct.prop <- colSums(counts.filter[chrY.genes,meta.filter$cell])/colSums(counts.filter[,meta.filter$cell])
temp.em.sex <- meta.filter.mod %>% mutate(sex=ifelse(chrY.RSEM.sum < exp.min, "F","Undef")) %>% mutate(sex=ifelse(chrY.RSEM.sum > exp.max, "M",sex)) %>% group_by(embryo,sex) %>% summarise(nC=n()) %>% group_by(embryo) %>% mutate(prop=nC/sum(nC),sum_nC=sum(nC)) %>% group_by(embryo) %>% top_n(1,prop)  %>% mutate(sex=ifelse(prop > exp.prop,sex,"Undef"))%>% unique()
meta.filter.mod <- meta.filter.mod %>% inner_join(temp.em.sex %>% select(embryo,sex),by="embryo")

saveRDS(meta.filter.mod,file=paste0("tmp_data/",TD,"/mouse_coseq_smt2.meta.filter.mod.withsex.rds"))

