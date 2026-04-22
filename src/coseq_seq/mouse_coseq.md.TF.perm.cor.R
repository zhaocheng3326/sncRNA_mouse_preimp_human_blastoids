#' ---
#' title: permutation for module miRNA - TF paris
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


#' coseq smt2 part
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

#' loading ME for coseq
md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
miRNA.MEs <-md.miRNA.list$MEs[,colnames(md.miRNA.list$MEs) %>% setdiff(c("grey"))]
miRNA.MEs.anno <- md.miRNA.list$md.anno
colnames(miRNA.MEs ) <- (miRNA.MEs.anno %>% tibble::column_to_rownames("miRNA_mdc"))[colnames(miRNA.MEs ),"miRNA_md"]


#' loading cor between md and TF and miRNA modules
md_miRNA_TF.cor.out <- readRDS(paste0("tmp_data/",TD,"/md_miRNA_TF.cor.out.rds"))
coseq.md.exp <- as.data.frame(t(miRNA.MEs))[,sel_cells]
coseq.TF.exp <- coseq.gene.exp[unique(md_miRNA_TF.cor.out$TF),sel_cells]



#' generate sampling (miRNA module order/ cell order/ TF order)
if (file.exists(paste0("tmp_data/",TD,"/coseq.md.TF.perm.set.rds"))) {
  perm_list <- readRDS(paste0("tmp_data/",TD,"/coseq.md.TF.perm.set.rds"))
}else{
  #' don't overwrite
  md_name <- rownames(coseq.md.exp )
  TF_name <- rownames(coseq.TF.exp)
  cell_name <- sel_cells
  perm_list <- list()
  for (n in 1:1000 ) {
    print(n)
    perm_list[[paste0("perm_",n)]] <- list()
    perm_list[[paste0("perm_",n)]][["md_name"]] <- sample(md_name)
    perm_list[[paste0("perm_",n)]][["TF_name"]] <- sample(TF_name)
    perm_list[[paste0("perm_",n)]][["cell_name"]] <- sample(cell_name)
    perm_list[[paste0("perm_",n)]][["tag"]] <- paste0("perm_",n)
  }
  saveRDS(perm_list,paste0("tmp_data/",TD,"/coseq.md.TF.perm.set.rds"))
}



if (file.exists(paste0("tmp_data/",TD,"/","md_miRNA_TF.perm.perf.list.rds"))) {
  md_miRNA_TF.perm.perf.list <- readRDS(paste0("tmp_data/",TD,"/","md_miRNA_TF.perm.perf.list.rds"))
}else{
  left_right <- md_miRNA_TF.cor.out %>% select(TF,miRNA_md,TF) %>% unique()%>% setNames(c("left","right")) 
  md_miRNA_TF.perm.perf.list <- list()
  md_miRNA_TF.perm.perf.list <- list()
  for (n in 1:1000) {
    print(n)
    md_miRNA_TF.perm.perf.list[[paste0("perm_",n)]] <- Fun_para_md_TF_cor(perm_list,coseq.TF.exp,coseq.md.exp,left_right,paste0("perm_",n))
  }
  saveRDS(md_miRNA_TF.perm.perf.list,paste0("tmp_data/",TD,"/","md_miRNA_TF.perm.perf.list.rds"))
}


#' check the sig ones
if (file.exists(paste0("tmp_data/",TD,"/coseq.md_miRNA.TF.cor.out.withPermutation.pvalue.rds"))) {
  md_miRNA_TF.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.md_miRNA.TF.cor.out.withPermutation.pvalue.rds"))
}else{
  md_miRNA_TF.cor.pv <- md_miRNA_TF.cor.out  %>%select(miRNA_md,TF,r) %>% mutate(nC=0) %>% mutate(r=round(r,5))
  
  for (n in names(md_miRNA_TF.perm.perf.list)) {
    print(n)
    md_miRNA_TF.cor.pv <- left_join(md_miRNA_TF.cor.out %>%select(miRNA_md,TF,r),md_miRNA_TF.perm.perf.list[[n]] %>% rename(perm_r=r,miRNA_md=md),by=c("miRNA_md","TF")) %>% mutate(perm_r=ifelse(is.na(perm_r),0,perm_r)) %>% filter((perm_r > r  & r > 0) | (perm_r < r & r < 0)) %>% mutate(nC_add=1) %>% select(-c(r,perm_r)) %>% right_join(md_miRNA_TF.cor.pv,by=c("miRNA_md","TF")) %>% mutate(nC_add=ifelse(is.na(nC_add),0,nC_add)) %>% mutate(nC=nC+nC_add) %>% select(-nC_add)
  }

  
  md_miRNA_TF.cor.pv.out <- md_miRNA_TF.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r > 0) %>% split(.,.$miRNA_md) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.) %>% bind_rows(md_miRNA_TF.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r < 0) %>% split(.,.$miRNA_md) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.))
  saveRDS(md_miRNA_TF.cor.pv.out,paste0("tmp_data/",TD,"/coseq.md_miRNA.TF.cor.out.withPermutation.pvalue.rds"))
}

md_miRNA_TF.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% group_by(TF) %>% summarise(nM=n_distinct(miRNA_md)) %>% arrange(desc(nM))
md_miRNA_TF.cor.pv.out %>% filter(r > 0 & FDR < 0.05) %>% group_by(TF) %>% summarise(nM=n_distinct(miRNA_md)) %>% arrange(desc(nM))

