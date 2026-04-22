#' ---
#' title: cor for gene and miRNA modules
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
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the full cell annotation
meta.filter <- meta.filter  %>% inner_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap %>% filter(EML!="unknown"),by="cell") %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))

#' loading ME for coseq
md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
md.gene.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_smt2",".md.list.rds"))

#' the rank for miRNA md
md.miRNA.list$kME.rank <- md.miRNA.list$modules %>% tbl_df() %>% rename(miRNA_name=gene_name) %>% gather(miRNA_mdc,kME,-c(miRNA_name,module,color)) %>% tbl_df() %>% mutate(miRNA_mdc=gsub("kME_","",miRNA_mdc)) %>% filter(miRNA_mdc==color) %>% select(miRNA_name,miRNA_mdc,kME) %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") %>% split(.,.$miRNA_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)
md.gene.list$kME.rank <- md.gene.list$modules %>% tbl_df()%>% gather(gene_mdc,kME,-c(gene_name,module,color)) %>% tbl_df() %>% mutate(gene_mdc=gsub("kME_","",gene_mdc)) %>% filter(gene_mdc==color) %>% select(gene_name,gene_mdc,kME) %>% inner_join(md.gene.list$md.anno,by="gene_mdc") %>% split(.,.$gene_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)


gene.MEs <- md.gene.list$MEs[,colnames(md.gene.list$MEs) %>% setdiff(c("grey"))]
miRNA.MEs <-md.miRNA.list$MEs[,colnames(md.miRNA.list$MEs) %>% setdiff(c("grey"))]

gene.MEs.anno <- md.gene.list$md.anno
miRNA.MEs.anno <- md.miRNA.list$md.anno

colnames(gene.MEs ) <- (gene.MEs.anno %>% tibble::column_to_rownames("gene_mdc"))[colnames(gene.MEs ),"gene_md"]
colnames(miRNA.MEs ) <- (miRNA.MEs.anno %>% tibble::column_to_rownames("miRNA_mdc"))[colnames(miRNA.MEs ),"miRNA_md"]

sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% pull(cell) %>% intersect(rownames(md.gene.list$MEs))%>% intersect(rownames(md.miRNA.list$MEs))

#' miRNA exp
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))[,sel_cells]
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))

#' module target genes
merged.target.pairs <- target.pairs %>% inner_join( md.miRNA.list$modules %>% select(gene_name,module) %>% tbl_df() %>% rename(miRNA=gene_name,miRNA_mdc=module),by="miRNA") %>% group_by(gene,miRNA_mdc) %>% summarise(nTargets=n_distinct(miRNA)) %>% ungroup()

if (file.exists(paste0("tmp_data/",TD,"/md_miRNA_md_gene.cor.out.rds"))) {
  md_miRNA_md_gene.cor.out <- readRDS(file=paste0("tmp_data/",TD,"/md_miRNA_md_gene.cor.out.rds"))

}else{
  #' calculate correlation between miRNA and gene modules
  left_data <- as.data.frame(t(miRNA.MEs[sel_cells,]))
  right_data <- as.data.frame(t(gene.MEs[sel_cells,]))
  left_right <- expand.grid(rownames(left_data),rownames(right_data)) %>% as.data.frame() %>% setNames(c("left","right"))
  md_miRNA_gene.cor <- FunPairedCor(left_data ,right_data,left_right,method = "spearman") %>% ungroup() %>% rename(miRNA_md=left,gene_md=right)
  
  md_miRNA_md_gene.cor.out <- md_miRNA_gene.cor %>% mutate(fdr=p.adjust(p)) %>% inner_join(miRNA.MEs.anno ,by="miRNA_md")%>% inner_join(gene.MEs.anno ,by="gene_md") %>% arrange(r)
  
  #' for blastocyst only stages (ICM vs TE)
  temp_sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% filter(EML %in% c("ICM","TE")) %>% pull(cell) 
  left_data <- as.data.frame(t(miRNA.MEs[temp_sel_cells,]))
  right_data <- as.data.frame(t(gene.MEs[temp_sel_cells,]))
  left_right <- expand.grid(rownames(left_data),rownames(right_data)) %>% as.data.frame() %>% setNames(c("left","right"))
  md_miRNA_gene.blast.cor <- FunPairedCor(left_data ,right_data,left_right,method = "spearman") %>% ungroup() %>% rename(miRNA_md=left,gene_md=right)
  md_miRNA_gene.blast.cor.out <- md_miRNA_gene.blast.cor %>% mutate(fdr=p.adjust(p)) %>% inner_join(miRNA.MEs.anno ,by="miRNA_md")%>% inner_join(gene.MEs.anno ,by="gene_md") %>% arrange(r)
  
  
  #' saving object
  saveRDS(md_miRNA_md_gene.cor.out ,file=paste0("tmp_data/",TD,"/md_miRNA_md_gene.cor.out.rds"))
  saveRDS(md_miRNA_gene.blast.cor.out ,file=paste0("tmp_data/",TD,"/md_miRNA_md_gene.blast.cor.out.rds"))
}

