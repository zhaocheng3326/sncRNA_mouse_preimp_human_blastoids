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

#' average expression
miRNA.ave.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.ave.exp.rds"))

temp <-  read.delim("big_doc/sc_smallRNA_annotation/Mouse/M.musculus.tsv",stringsAsFactors = F,head=T,sep="\t") %>% tbl_df()
TF_miRNA <- readRDS(paste0("tmp_data/",TD,"/TF_miRNA.rds"))

#' check TF
md_miRNA_TF.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.md_miRNA.TF.cor.out.withPermutation.pvalue.rds"))
md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
md.gene.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_smt2",".md.list.rds"))




TF_miRNA <- readRDS(paste0("tmp_data/",TD,"/TF_miRNA.rds")) ## TF regulating miRNA
TF.exp <-  readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))[unique(TF_miRNA$TF),]
miRNA.bed <- readRDS(paste0("tmp_data/",TD,"/miRNA.bed.rds"))

md.miRNA.list$kME.rank <- md.miRNA.list$modules %>% tbl_df() %>% rename(miRNA_name=gene_name) %>% gather(miRNA_mdc,kME,-c(miRNA_name,module,color)) %>% tbl_df() %>% mutate(miRNA_mdc=gsub("kME_","",miRNA_mdc)) %>% filter(miRNA_mdc==color) %>% select(miRNA_name,miRNA_mdc,kME) %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") %>% split(.,.$miRNA_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME") %>% mutate(quant_rank_kME=ntile(-1*kME,5)) }) %>% do.call("bind_rows",.)
md.gene.list$kME.rank <- md.gene.list$modules %>% tbl_df()%>% gather(gene_mdc,kME,-c(gene_name,module,color)) %>% tbl_df() %>% mutate(gene_mdc=gsub("kME_","",gene_mdc)) %>% filter(gene_mdc==color) %>% select(gene_name,gene_mdc,kME) %>% inner_join(md.gene.list$md.anno,by="gene_mdc") %>% split(.,.$gene_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME") %>% mutate(quant_rank_kME=ntile(-1*kME,5)) }) %>% do.call("bind_rows",.)


miRNA_TF.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.TF.cor.out.withPermutation.pvalue.rds"))
miRNA.neg.sc.TF.module.enriched.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TF.module.enriched.rds"))
miRNA.pos.sc.TF.module.enriched.out <- readRDS(file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.pos.TF.module.enriched.rds"))

md.miRNA.TF.enrich.out <- miRNA.pos.sc.TF.module.enriched.out %>% bind_rows(miRNA.neg.sc.TF.module.enriched.out )


#' select miRNA modules and TF
md.miRNA.TF.enrich.out.detail<- md.miRNA.TF.enrich.out %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") %>% inner_join(md_miRNA_TF.cor.pv.out %>% mutate(TF_md_cor=r,TF_md_pvalue=pvalue,TF_md_fdr=FDR)  %>% select(miRNA_md,TF,TF_md_cor,TF_md_pvalue,TF_md_fdr),by=c("miRNA_md","TF")) 


md.miRNA.TF.enrich.out.filter <- md.miRNA.TF.enrich.out.detail %>% filter(TF_md_fdr < 0.05 & p_val_adj < 0.05)



temp <- md.miRNA.TF.enrich.out.filter %>% inner_join(md.miRNA.list$modules %>% rename(miRNA=gene_name,miRNA_mdc=module) %>% select(miRNA,miRNA_mdc) %>% tbl_df()%>% unique(),by=c("miRNA_mdc"),relationship = "many-to-many") %>% inner_join(miRNA_TF.cor.pv.out %>% filter(FDR < 0.05) %>% mutate(TF_miRNA_cor=r,TF_miRNA_FDR=FDR) %>% select(TF,miRNA,TF_miRNA_cor,TF_miRNA_FDR) %>% unique(),by=c("TF","miRNA")) 

temp.id <- temp %>% filter(miRNA_mdc=="blue" & TF=="Yy1" & TF_miRNA_FDR < 0.05) %>% select(miRNA,TF_miRNA_cor) %>% rename(V7=miRNA) %>% inner_join(miRNA.bed ,by="V7") %>% arrange(TF_miRNA_cor) %>% mutate(pos=paste0(V1,":",V2,"-",V3)) %>% pull(V7) %>%unique()
(miRNA.ave.exp %>% filter(gene %in% temp.id) %>% tibble::column_to_rownames("gene"))[,EML.od] %>% pheatmap(scale="row",cluster_cols = F)


(miRNA.ave.exp %>% filter(gene %in% temp.id) %>% tibble::column_to_rownames("gene"))[,EML.od] %>% pheatmap(scale="row",cluster_cols = F)
mature_prec.anno <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/small_prec.cor.txt",sep="\t",stringsAsFactors=F,header = F) %>% tbl_df()

#' check 
####intersectBed -a Yy1_8C.merged_peaks.narrowPeak -b ~/My_project/mouse_smallseq_preimp/big_doc/sc_smallRNA_annotation/Mouse/prec.miRNA.up5k.bed  -wa -wb  -nonamecheck > a
#read.delim("/home/chenzh/My_project/mouse_smallseq_preimp/tmp_data/Yy1_CUTRUN_8C_GSM7122977/macs_Yy1_8C_merged/a",head=F) %>% mutate(V2=V17) %>% select(V2)%>% unique() %>% inner_join(mature_prec.anno,by="V2") %>% tbl_df() %>% pull(V1) %>% intersect(temp.id)

#intersectBed -a Yy1_ES_peaks.narrowPeak  -b ~/My_project/mouse_smallseq_preimp/big_doc/sc_smallRNA_annotation/Mouse/prec.miRNA.up5k.bed  -wa -wb  -nonamecheck > a
#read.delim("~/My_project/mouse_smallseq_preimp/tmp_data/Yy1_CHIP_GSE99518_CELL/a",head=F) %>% mutate(V2=V17) %>% select(V2)%>% unique() %>% inner_join(mature_prec.anno,by="V2") %>% tbl_df() %>% pull(V1) %>% intersect(temp.id)



