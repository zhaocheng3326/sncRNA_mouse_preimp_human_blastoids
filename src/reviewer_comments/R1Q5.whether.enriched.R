#' ---
#' title: answer Dlk1-Dio3 enriched or not
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

###The Dlk1-Dio3 region encodes an unusually large number of miRNA species, so counting species alone can artificially inflate apparent enrichment relative to other loci. Without a null model (e.g., permutation preserving family structure), the 79.5% value reported in Figure 3a,c may reflect cluster size rather than true ICM specificity. Support this claim beyond descriptive visualization requires a statistical enrichment test and abundance-containing metrics.

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
  library(glmnet)
  #library(pROC)
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


para_cal=TRUE
if (para_cal) {
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  numCores <- 11
  registerDoParallel(numCores)
}


#' C2MC12MC
C2C12MC.ov <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/C2MC_C12MC_ov_small_anno.bed",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V4 %in% c("C2MC","C12MC")) %>% filter(V12=="miRNA") %>% select(V4,V11,V12) %>% unique()
C12MC.miRNA.ID <- C2C12MC.ov %>% filter(V4=="C12MC",V12=="miRNA") %>% pull(V11) %>% unique()
C2MC.miRNA.ID <- C2C12MC.ov %>% filter(V4=="C2MC",V12=="miRNA") %>% pull(V11) %>% unique()


#loading marker 
miRNA.mk.out <- readRDS(file=paste0("tmp_data/",TD,"/","allCells.miRNA.mk.rds"))
ICM.mk.miRNA <- miRNA.mk.out %>% filter(cluster=="ICM") %>% pull(gene) %>% unique()
ICM_TE.mk.miRNA <- miRNA.mk.out %>% filter(cluster=="ICM_TE") %>% pull(gene) %>% unique()

#' expressed miRNA
miRNA.expressed <- rownames(readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds")))

#' check the C12MC cluster
temp <- data.frame(TSet=length(intersect(ICM.mk.miRNA,C12MC.miRNA.ID)), BgSet=length(ICM.mk.miRNA),TBg=length(intersect(miRNA.expressed,C12MC.miRNA.ID)),Bg=length(miRNA.expressed)) 
temp.fish.cal <-apply(temp,1,function(x){fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE),alternative = "great")} %>% return())
temp.fish.cal[[1]]$p.value
temp.fish.cal[[1]]$estimate


#' check the C2MC cluster
temp <- data.frame(TSet=length(intersect(ICM_TE.mk.miRNA,C2MC.miRNA.ID)), BgSet=length(ICM_TE.mk.miRNA),TBg=length(intersect(miRNA.expressed,C2MC.miRNA.ID)),Bg=length(miRNA.expressed)) 
temp.fish.cal <-apply(temp,1,function(x){fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE),alternative = "great")} %>% return())
temp.fish.cal[[1]]$p.value
temp.fish.cal[[1]]$estimate




#‘ lineage marker permutation
#' generate sample meta.filter
if (file.exists(paste0("tmp_data/",TD,"/","R1Q5.perm.ICM_TE.marker.rds"))) {
  perm.mk.list <- readRDS(paste0("tmp_data/",TD,"/","R1Q5.perm.ICM_TE.marker.rds"))
}else{
  # must save to repeat
  perm.mk.list <- list()
  mk.type <- data.frame(miRNA=miRNA.expressed,cluster="other") %>% tbl_df()%>% rows_update(miRNA.mk.out %>% rename(miRNA=gene) %>% select(miRNA,cluster),by="miRNA")
  for (n in 1:1000) {
    temp.mk.type <- mk.type %>% mutate(perm=paste0("perm_",n))
    temp.mk.type$cluster=sample(temp.mk.type$cluster)
    perm.mk.list[[paste0("perm_",n)]] <- temp.mk.type
  }
  saveRDS(perm.mk.list ,paste0("tmp_data/",TD,"/","R1Q5.perm.ICM_TE.marker.rds"))
}
#' get the permutation p values
temp.real.C12MC.number <- length(intersect(ICM.mk.miRNA,C12MC.miRNA.ID))
perm.mk.list %>% do.call("bind_rows",.) %>% filter(cluster=="ICM") %>% filter(miRNA %in% C12MC.miRNA.ID ) %>% group_by(perm) %>% summarise(nMiRNA=n_distinct(miRNA)) %>% filter(nMiRNA >=temp.real.C12MC.number  )

temp.real.C2MC.number <- length(intersect(ICM_TE.mk.miRNA,C2MC.miRNA.ID))
perm.mk.list %>% do.call("bind_rows",.) %>% filter(cluster=="ICM_TE") %>% filter(miRNA %in% temp.real.C2MC.number  ) %>% group_by(perm) %>% summarise(nMiRNA=n_distinct(miRNA)) %>% filter(nMiRNA >=temp.real.C2MC.number  )

