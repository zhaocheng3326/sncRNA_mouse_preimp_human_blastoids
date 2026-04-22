#' ---
#' title: permutation for miRNA - gene paris (only blastocyst cells)
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
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))[,sel_cells] %>% as.data.frame()
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))[,sel_cells]%>% as.data.frame()

#' generate sampling (miRNA order/ cell order/ gene order)
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.gene.perm.set.rds"))) {
  perm_list <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.gene.perm.set.rds"))
}

rerun=FALSE
if (rerun) {
  left_right <- target.pairs%>% select(miRNA,gene) %>% setNames(c("left","right"))
  x=1:1000;
  for (ns in split(x, ceiling(seq_along(x)/100))) {
    nmin <- min(ns)
    nmax <- max(ns)
    OB.miRNA.gene.perm.perf.list <- list()
    for (n in nmin:nmax) {
      print(n)
      OB.miRNA.gene.perm.perf.list[[paste0("perm_",n)]] <- Fun_para_miRNA_gene_cor(perm_list,coseq.miRNA.exp,coseq.gene.exp,left_right,paste0("perm_",n),blast_sel_cells )
    } ###  ### too big, need to split
    saveRDS(OB.miRNA.gene.perm.perf.list,paste0("tmp_data/",TD,"/","temp.perm.OB.P",nmin,"_to_",nmax,".rds"))
  }
}


#' check the sig ones
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.withPermutation.pvalue.rds"))) {
  miRNA_gene.blast.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.withPermutation.pvalue.rds"))
}else{
  miRNA_gene.blast.cor <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.rds"))
  x=1:1000;
  miRNA_gene.blast.cor.pv <- miRNA_gene.blast.cor %>%select(miRNA,gene,r) %>% mutate(nC=0) %>% mutate(r=round(r,5))
  for (ns in split(x, ceiling(seq_along(x)/100))) {
    nmin <- min(ns)
    nmax <- max(ns)
    
    miRNA.gene.perm.perf.list <- readRDS(paste0("tmp_data/",TD,"/","temp.perm.OB.P",nmin,"_to_",nmax,".rds"))
    for (n in names(miRNA.gene.perm.perf.list)) {
      print(n)
      miRNA_gene.blast.cor.pv <- left_join(miRNA_gene.blast.cor %>%select(miRNA,gene,r),miRNA.gene.perm.perf.list[[n]] %>% rename(perm_r=r),by=c("miRNA","gene")) %>% mutate(perm_r=ifelse(is.na(perm_r),0,perm_r)) %>% filter((perm_r > r  & r > 0) | (perm_r < r & r < 0)) %>% mutate(nC_add=1) %>% select(-c(r,perm_r)) %>% right_join(miRNA_gene.blast.cor.pv,by=c("miRNA","gene")) %>% mutate(nC_add=ifelse(is.na(nC_add),0,nC_add)) %>% mutate(nC=nC+nC_add) %>% select(-nC_add)
    }
  }
  
  miRNA_gene.blast.cor.pv.out <- miRNA_gene.blast.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r > 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.) %>% bind_rows(miRNA_gene.blast.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r < 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.))
  saveRDS(miRNA_gene.blast.cor.pv.out,paste0("tmp_data/",TD,"/coseq.miRNA.gene.blast.cor.out.withPermutation.pvalue.rds"))
}