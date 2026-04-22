#' ---
#' title: feature of miRNA modules
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

TF_miRNA <- readRDS(paste0("tmp_data/",TD,"/TF_miRNA.rds"))
miRNA.bed <- readRDS(paste0("tmp_data/",TD,"/miRNA.bed.rds"))
# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=recode(EML,"prelineage"="L8CM")) %>% select(cell,RNA_EML,small_EML,EML)  #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the full cell annotation
meta.filter <- meta.filter  %>% inner_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap %>% filter(EML!="unknown"),by="cell")%>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))

#' coseq smt2 part
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

#' loading ME for coseq
md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
miRNA.MEs <-md.miRNA.list$MEs[,colnames(md.miRNA.list$MEs) %>% setdiff(c("grey"))]
miRNA.MEs.anno <- md.miRNA.list$md.anno
colnames(miRNA.MEs ) <- (miRNA.MEs.anno %>% tibble::column_to_rownames("miRNA_mdc"))[colnames(miRNA.MEs ),"miRNA_md"]

#' using the cells have the same RNA_EML and small_EML
sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% pull(cell) %>% intersect(readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds"))$cell)%>% intersect(rownames(md.miRNA.list$MEs))
coseq.miRNA.exp <- coseq.miRNA.exp[,sel_cells]


#md.miRNA.list$modules %>% gather(miRNA_mdc,kME,-c(gene_name,module,color)) %>% tbl_df() %>% mutate(miRNA_mdc=gsub("kME_","",miRNA_mdc)) %>% filter(miRNA_mdc==color)
#' TF enrichment in miRNA modules
md.miRNA.TF.enrich <- list()

for (m in unique(md.miRNA.list$modules$module)) {
  md.miRNA.TF.enrich[[m]] <- md.miRNA.list$modules %>% filter(module==m) %>% pull(gene_name) %>% FunMiRNA_TF_fisher(TF_miRNA) %>% mutate(miRNA_mdc=m)
}

md.miRNA.TF.enrich.out <- md.miRNA.TF.enrich %>% do.call("bind_rows",.) %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") #%>% mutate(enrichment=ifelse(p_val_adj < 0.05 & od  >=1,"Sig","notSig"))
md.miRNA.TF.enrich.out %>% filter(pvalue < 0.05) %>% arrange(p_val_adj)


#' TF correlation with miRNA modules
left_data <-  as.data.frame(t(miRNA.MEs[sel_cells,]))
right_data <- coseq.gene.exp[,sel_cells]
left_right <- md.miRNA.TF.enrich.out  %>% arrange(p_val_adj) %>% select(miRNA_md,TF) %>% unique()#%>% filter(enrichment=="Sig")

md_miRNA_TF.cor <- FunPairedCor(left_data ,right_data,left_right,method = "spearman") %>% ungroup() %>% rename(miRNA_md=left,TF=right)
md_miRNA_TF.cor.out <- md_miRNA_TF.cor %>% mutate(fdr=p.adjust(p)) %>% inner_join(miRNA.MEs.anno ,by="miRNA_md")%>% arrange(desc(abs(r)))

md_miRNA_TF.cor.out %>% filter(fdr < 0.05 ) %>% arrange(desc(abs(r))) 
md.miRNA.TF.enrich.out.filter <- md.miRNA.TF.enrich.out %>% inner_join(md_miRNA_TF.cor.out %>% mutate(TF_cor=r,TF_pvalue=p,TF_fdr=fdr)  %>% select(miRNA_md,miRNA_mdc,TF,TF_cor,TF_pvalue,TF_fdr),by=c("miRNA_md","miRNA_mdc","TF"))%>% filter(TF_pvalue < 0.05 & pvalue < 0.05)

#' TF correlation with miRNAs
left_data <-  coseq.miRNA.exp[,sel_cells] %>% as.data.frame()
right_data <- coseq.gene.exp[,sel_cells]
left_right <- TF_miRNA %>% select(miRNA,TF) %>% unique()

miRNA_TF.cor <- FunPairedCor(left_data ,right_data,left_right,method = "spearman") %>% ungroup() %>% rename(miRNA=left,TF=right)

rewrite=FALSE
#' saving object
if (rewrite) {
  saveRDS(md_miRNA_TF.cor.out,file=paste0("tmp_data/",TD,"/md_miRNA_TF.cor.out.rds"))
  saveRDS(md.miRNA.TF.enrich.out ,file=paste0("tmp_data/",TD,"/md.miRNA.TF.enrich.rds"))
  saveRDS(miRNA_TF.cor,file=paste0("tmp_data/",TD,"/miRNA.TF.cor.out.rds"))
}
