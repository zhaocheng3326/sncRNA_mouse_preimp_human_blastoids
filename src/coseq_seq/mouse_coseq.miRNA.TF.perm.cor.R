#' ---
#' title: permutation for miRNA - TF paris
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
  left_right <- TF_miRNA%>% select(miRNA,TF)  %>% unique()%>% setNames(c("left","right")) 
  x=1:1000;
  for (ns in split(x, ceiling(seq_along(x)/100))) {
    nmin <- min(ns)
    nmax <- max(ns)
    miRNA.TF.perm.perf.list <- list()
    for (n in nmin:nmax) {
      print(n)
      miRNA.TF.perm.perf.list[[paste0("perm_",n)]] <- Fun_para_miRNA_gene_cor(perm_list,coseq.miRNA.exp,coseq.gene.exp,left_right,paste0("perm_",n))
    } ###  ### too big, need to split
    saveRDS(miRNA.TF.perm.perf.list,paste0("tmp_data/",TD,"/","temp.perm.TF.miRNA.P",nmin,"_to_",nmax,".rds"))
  }
}


#' check the sig ones
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.TF.cor.out.withPermutation.pvalue.rds"))) {
  miRNA_TF.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.TF.cor.out.withPermutation.pvalue.rds"))
}else{
  miRNA_TF.cor <- readRDS(paste0("tmp_data/",TD,"/miRNA.TF.cor.out.rds"))
  x=1:1000;
  miRNA_TF.cor.pv <- miRNA_TF.cor %>%select(miRNA,TF,r) %>% mutate(nC=0) %>% mutate(r=round(r,5))
  for (ns in split(x, ceiling(seq_along(x)/100))) {
    nmin <- min(ns)
    nmax <- max(ns)
    
    miRNA.TF.perm.perf.list <- readRDS(paste0("tmp_data/",TD,"/","temp.perm.TF.miRNA.P",nmin,"_to_",nmax,".rds"))
    for (n in names(miRNA.TF.perm.perf.list)) {
      print(n)
      miRNA_TF.cor.pv <- left_join(miRNA_TF.cor %>%select(miRNA,TF,r),miRNA.TF.perm.perf.list[[n]] %>% rename(perm_r=r,TF=gene),by=c("miRNA","TF")) %>% mutate(perm_r=ifelse(is.na(perm_r),0,perm_r)) %>% filter((perm_r > r  & r > 0) | (perm_r < r & r < 0)) %>% mutate(nC_add=1) %>% select(-c(r,perm_r)) %>% right_join(miRNA_TF.cor.pv,by=c("miRNA","TF")) %>% mutate(nC_add=ifelse(is.na(nC_add),0,nC_add)) %>% mutate(nC=nC+nC_add) %>% select(-nC_add)
    }
  }
  
  miRNA_TF.cor.pv.out <- miRNA_TF.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r > 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.) %>% bind_rows(miRNA_TF.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r < 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.))
  saveRDS(miRNA_TF.cor.pv.out,paste0("tmp_data/",TD,"/coseq.miRNA.TF.cor.out.withPermutation.pvalue.rds"))
}

miRNA_TF.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% group_by(TF) %>% summarise(nM=n_distinct(miRNA)) %>% arrange(desc(nM))
miRNA_TF.cor.pv.out %>% filter(r > 0 & FDR < 0.05) %>% group_by(TF) %>% summarise(nM=n_distinct(miRNA)) %>% arrange(desc(nM))


#' check the enrichment of targeted genes on each gene module
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TG.pt.rds"))) {
  miRNA.neg.sc.TF.module.enriched.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TF.module.enriched.rds"))
  miRNA.pos.sc.TF.module.enriched.out <- readRDS(file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.pos.TF.module.enriched.rds"))
}else{
  md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
  
  #' the rank for miRNA md
  md.miRNA.list$kME.rank <- md.miRNA.list$modules %>% tbl_df() %>% rename(miRNA_name=gene_name) %>% gather(miRNA_mdc,kME,-c(miRNA_name,module,color)) %>% tbl_df() %>% mutate(miRNA_mdc=gsub("kME_","",miRNA_mdc)) %>% filter(miRNA_mdc==color) %>% select(miRNA_name,miRNA_mdc,kME) %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") %>% split(.,.$miRNA_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)
  
 
  #' neg sig cor ones
  temp.neg.sig.target.pairs <- miRNA_TF.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% select(miRNA,TF) %>% unique()
  temp.list <- list()
  for (m in unique(md.miRNA.list$modules$module)) {
    print(m)
    temp.list[[m]] <- md.miRNA.list$kME.rank %>% filter(miRNA_mdc==m)  %>% pull(miRNA_name) %>% FunMiRNA_TF_fisher(temp.neg.sig.target.pairs) %>% mutate(miRNA_mdc=m)
  }
  miRNA.neg.sc.TF.module.enriched.out <- temp.list %>% do.call("bind_rows",.) %>% arrange(p_val_adj)
  
  #' pos sig cor ones
  temp.pos.sig.target.pairs <- miRNA_TF.cor.pv.out %>% filter(r > 0 & FDR < 0.05) %>% select(miRNA,TF) %>% unique()
  temp.list <- list()
  for (m in unique(md.miRNA.list$modules$module)) {
    print(m)
    temp.list[[m]] <- md.miRNA.list$kME.rank %>% filter(miRNA_mdc==m)  %>% pull(miRNA_name) %>% FunMiRNA_TF_fisher(temp.pos.sig.target.pairs ) %>% mutate(miRNA_mdc=m)
  }
  miRNA.pos.sc.TF.module.enriched.out <- temp.list %>% do.call("bind_rows",.) %>% arrange(p_val_adj)
  
  
  saveRDS(miRNA.neg.sc.TF.module.enriched.out,file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TF.module.enriched.rds"))
  saveRDS(miRNA.pos.sc.TF.module.enriched.out,file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.pos.TF.module.enriched.rds"))
}





