#' ---
#' title: feature/ pathway/GO enrichment analysis for module genes
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
  library(clusterProfiler)
  library(topGO)
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


#' loading ME for coseq
md.gene.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_smt2",".md.list.rds"))

#' the rank for miRNA md
md.gene.list$kME.rank <- md.gene.list$modules %>% tbl_df()%>% gather(gene_mdc,kME,-c(gene_name,module,color)) %>% tbl_df() %>% mutate(gene_mdc=gsub("kME_","",gene_mdc)) %>% filter(gene_mdc==color) %>% select(gene_name,gene_mdc,kME) %>% inner_join(md.gene.list$md.anno,by="gene_mdc") %>% split(.,.$gene_mdc) %>% lapply(function(x){x %>% arrange(desc(kME)) %>% tibble::rowid_to_column("rank_kME")}) %>% do.call("bind_rows",.)


#' pathway annotaiton
ALL_gene <- rownames( readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.counts.filter.rds")))
geneID2GO=inverseList(readMappings(file = "big_doc/mouse.GO2geneID_ALL.refseq.map"))
GS_db =  readRDS("big_doc/Mus.msigdbr.rds")


gene.md.pt.results <- list()

if (file.exists(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.md.pt.results.rds"))) {
  gene.md.pt.results <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.md.pt.results.rds"))
}else{
  for (g in unique(md.gene.list$modules$module)) {
    print(g)
    gene.md.pt.results[[g]] <- list()
    temp.out <- list()
    
    temp.genes <- md.gene.list$modules %>% filter(module==g) %>% pull(gene_name) %>% unique()
    temp.out[["gene"]] <- temp.genes
    #' GO enrichment analysis
    temp.out[["GO"]] <- topGO_enrichment2(temp.genes,ALL_gene,0.05,geneID2GO)
    
    for (tp in c("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS")) { #"CGP","CP", "TFT:GTRD"
      temp.path <- GS_db %>% filter(gs_subcat == tp) %>% mutate(TERM=gs_name,GENE=gene_symbol) %>% select(TERM,GENE) %>% unique() %>% as.data.frame()
      suppressMessages(temp.out[[tp]]  <- enricher( temp.genes, pvalueCutoff = 1, pAdjustMethod = "BH",ALL_gene ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
      
    }
    #' enrichR 2024 wiki
    temp.path <- read.delim("~/Genome_new/enrichr_library/modified/WikiPathways_2024_Mouse.mod.txt",stringsAsFactors = F,head=F)%>% setNames(c("TERM","GENE")) %>% unique() %>% as.data.frame()
    suppressMessages(temp.out$enrichR_wiki  <- enricher( temp.genes, pvalueCutoff = 1, pAdjustMethod = "BH",ALL_gene ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
    
    gene.md.pt.results[[g]] <- temp.out
  }
  #' saving object
  saveRDS(gene.md.pt.results,paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.md.pt.results.rds"))
}


