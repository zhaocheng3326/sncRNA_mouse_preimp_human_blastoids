#' ---
#' title: permutation for miRNA - gene paris
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
}else{
  #' don't overwrite
  miRNA_name <- rownames(coseq.miRNA.exp)
  gene_name <- rownames(coseq.gene.exp)
  cell_name <- sel_cells
  perm_list <- list()
  for (n in 1:1000 ) {
    print(n)
    perm_list[[paste0("perm_",n)]] <- list()
    perm_list[[paste0("perm_",n)]][["miRNA_name"]] <- sample(miRNA_name)
    perm_list[[paste0("perm_",n)]][["gene_name"]] <- sample(gene_name)
    perm_list[[paste0("perm_",n)]][["cell_name"]] <- sample(cell_name)
    perm_list[[paste0("perm_",n)]][["tag"]] <- paste0("perm_",n)
  }
  saveRDS(perm_list,paste0("tmp_data/",TD,"/coseq.miRNA.gene.perm.set.rds"))
}


#' check the sig ones
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.withPermutation.pvalue.rds"))) {
  miRNA_gene.cor.pv.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.withPermutation.pvalue.rds"))
}else{
  miRNA_gene.cor <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.rds"))
  x=1:1000;
  miRNA_gene.cor.pv <- miRNA_gene.cor %>%select(miRNA,gene,r) %>% mutate(nC=0) %>% mutate(r=round(r,5))
  for (ns in split(x, ceiling(seq_along(x)/100))) {
    nmin <- min(ns)
    nmax <- max(ns)
    
    miRNA.gene.perm.perf.list <- readRDS(paste0("tmp_data/",TD,"/","temp.perm.P",nmin,"_to_",nmax,".rds"))
    for (n in names(miRNA.gene.perm.perf.list)) {
      print(n)
      miRNA_gene.cor.pv <- left_join(miRNA_gene.cor %>%select(miRNA,gene,r),miRNA.gene.perm.perf.list[[n]] %>% rename(perm_r=r),by=c("miRNA","gene")) %>% mutate(perm_r=ifelse(is.na(perm_r),0,perm_r)) %>% filter((perm_r > r  & r > 0) | (perm_r < r & r < 0)) %>% mutate(nC_add=1) %>% select(-c(r,perm_r)) %>% right_join(miRNA_gene.cor.pv,by=c("miRNA","gene")) %>% mutate(nC_add=ifelse(is.na(nC_add),0,nC_add)) %>% mutate(nC=nC+nC_add) %>% select(-nC_add)
    }
  }
  
  miRNA_gene.cor.pv.out <- miRNA_gene.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r > 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.) %>% bind_rows(miRNA_gene.cor.pv %>% mutate(pvalue=nC/1000) %>% filter(r < 0) %>% split(.,.$miRNA) %>% lapply(function(x) {x %>% mutate(FDR=p.adjust(pvalue,method="fdr"))}) %>% do.call("bind_rows",.))
  saveRDS(miRNA_gene.cor.pv.out,paste0("tmp_data/",TD,"/coseq.miRNA.gene.cor.out.withPermutation.pvalue.rds"))
}

miRNA_gene.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% group_by(miRNA) %>% summarise(nTG=n_distinct(gene)) %>% arrange(desc(nTG))
#' check the overlap of let-7 family
print(
  miRNA_gene.cor.pv.out %>% filter(FDR < 0.05 & r < 0) %>%  filter(miRNA %in% c("mmu-let-7e-5p","mmu-let-7c-5p","mmu-let-7a-5p","mmu-let-7d-5p")) %>% mutate(miRNA=gsub("mmu-","",miRNA)) %>% split(.,.$miRNA) %>% lapply(function(x){x$gene}) %>% ggvenn::ggvenn(text_size=3,set_name_size=4)+ggtitle("Sig neg-correlated target genes\n(permutation, r<0, fdr < 0.05)")+FunTitle()
)

#' check the enrichment of targeted genes on each gene module
if (file.exists(paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TG.pt.rds"))) {
  miRNA.neg.sc.tg.module.enriched.out <- readRDS(paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TG.module.enriched.rds"))
  miRNA.pos.sc.tg.module.enriched.out <- readRDS(file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.pos.TG.module.enriched.rds"))
}else{
  md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
  md.gene.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_smt2",".md.list.rds"))
  
  #' the rank for miRNA md
  md.miRNA.list$kME.rank <- md.miRNA.list$modules %>% tbl_df() %>% rename(miRNA_name=gene_name) %>% gather(miRNA_mdc,kME,-c(miRNA_name,module,color)) %>% tbl_df() %>% mutate(miRNA_mdc=gsub("kME_","",miRNA_mdc)) %>% filter(miRNA_mdc==color) %>% select(miRNA_name,miRNA_mdc,kME) %>% inner_join(md.miRNA.list$md.anno,by="miRNA_mdc") %>% split(.,.$miRNA_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)
  
  md.gene.list$kME.rank <- md.gene.list$modules %>% tbl_df()%>% gather(gene_mdc,kME,-c(gene_name,module,color)) %>% tbl_df() %>% mutate(gene_mdc=gsub("kME_","",gene_mdc)) %>% filter(gene_mdc==color) %>% select(gene_name,gene_mdc,kME) %>% inner_join(md.gene.list$md.anno,by="gene_mdc") %>% split(.,.$gene_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)
  
  #' neg sig cor ones
  temp.neg.sig.target.pairs <- miRNA_gene.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% select(miRNA,gene) %>% unique()
  temp.list <- list()
  for (g in unique(md.gene.list$modules$module)) {
    print(g)
    temp.list[[g]] <- md.gene.list$kME.rank %>% filter(gene_mdc==g)  %>% pull(gene_name) %>% FunGeneSet_miRNA_fisher(temp.neg.sig.target.pairs) %>% mutate(gene_mdc=g)
  }
  miRNA.neg.sc.tg.module.enriched.out <- temp.list %>% do.call("bind_rows",.) %>% left_join(md.miRNA.list$kME.rank %>% rename(miRNA=miRNA_name),by="miRNA") %>% arrange(p_val_adj)
  
  #' pos sig cor ones
  temp.pos.sig.target.pairs <- miRNA_gene.cor.pv.out %>% filter(r > 0 & FDR < 0.05) %>% select(miRNA,gene) %>% unique()
  temp.list <- list()
  for (g in unique(md.gene.list$modules$module)) {
    print(g)
    temp.list[[g]] <- md.gene.list$kME.rank %>% filter(gene_mdc==g)  %>% pull(gene_name) %>% FunGeneSet_miRNA_fisher( temp.pos.sig.target.pairs ) %>% mutate(gene_mdc=g)
  }
  miRNA.pos.sc.tg.module.enriched.out <- temp.list %>% do.call("bind_rows",.) %>% left_join(md.miRNA.list$kME.rank %>% rename(miRNA=miRNA_name),by="miRNA") %>% arrange(p_val_adj)
  
 saveRDS(miRNA.neg.sc.tg.module.enriched.out,file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.neg.TG.module.enriched.rds"))
 saveRDS(miRNA.pos.sc.tg.module.enriched.out,file=paste0("tmp_data/",TD,"/coseq.miRNA.sig.pos.TG.module.enriched.rds"))
}


miRNA_gene.cor.pv.out %>% filter(r < 0 & FDR < 0.05) %>% inner_join( md.miRNA.list$modules %>% select(gene_name,module) %>% rename(miRNA=gene_name,miRNA_mdc=module),by="miRNA") %>% inner_join( md.gene.list$modules %>% select(gene_name,module) %>% rename(gene=gene_name,gene_mdc=module),by="gene") %>% group_by(miRNA) %>% summarise(nTG=n_distinct(gene),nGene_mdc=n_distinct(gene_mdc)) %>% arrange(desc(nTG)) %>% mutate(ratio=nTG/nGene_mdc)
miRNA_gene.cor.pv.out %>% filter(FDR < 0.05 & r<0) %>% group_by(miRNA) %>% summarise(nTG=n_distinct(gene)) %>% arrange(desc(nTG))
