#' ---
#' title: "coseq raw dataset human blasotids"
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


#' loading previous annotation of JP
counts <- list()
metas <- list()



# smt2 part for split seq
pj="HB_co_smt2"

temp <- read.delim("doc/HB_ME_July_2025/HB_ME_July_2025_smt2.txt",sep="\t",stringsAsFactors = F) %>% tbl_df() %>% tbl_df() %>% select(Cell..:Stage) %>% select(-Batch)%>% setNames(c("cell","embryo","Strain","stage")) %>% filter(Strain=="H9") %>% mutate(batch="HB_coseq",cell=paste0("Sample_",cell,"_0"),embryo=paste0("blastoids_",embryo),devTime=recode(stage,"16 cell"="16C","32 cell blast"="32C","32 cell morula"="32C","64 cell"="64C","8 cell (comp)"="8C","8 cell (uncomp)"="8C","2 cell"="2C","4 cell"="4C")) %>% select(-Strain)  %>% inner_join(read.delim("tmp_data/HB_coseq_smt2/multiqc_data/multiqc_star.txt",,stringsAsFactors = F) %>% tbl_df() %>% mutate(cell=Sample,mapped_percent=uniquely_mapped_percent+multimapped_percent) %>% select(cell,total_reads,mapped_percent,uniquely_mapped_percent),by="cell") 

counts[[pj]] <- read.delim("tmp_data/HB_coseq_smt2/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()           

metas[[pj]] <- temp %>% tbl_df() %>% mutate(seqType="smt2",cellType="EM")

counts[[pj]] <- counts[[pj]][,c("Gene",metas[[pj]]$cell)]
metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)                             




# Do the QC based on the MT.percent and nGene
gtf.anno <- read.delim( "/home/chenzh/Genome_new/Human/RefSeq/clean_chr_refdata/genes/gene.gtf.anno",stringsAsFactors = F,row.names = 5,head=F)
mt.gene <- gtf.anno$V6[gtf.anno$V1=="chrMT"]
ribo.gene <- c(gtf.anno$V6[grepl("^RPS",gtf.anno$V6)],gtf.anno$V6[grepl("^RPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPS",gtf.anno$V6)])


dup.gene <- names(table(gtf.anno$V6))[table(gtf.anno$V6 ) > 1]
GI <- unique(gtf.anno$V6)
for (n in c("HB_co_smt2")) {
  counts[[n]]$Gene <- gtf.anno[counts[[n]]$Gene,]$V6
  counts[[n]] <- (counts[[n]] %>% filter(Gene %in% dup.gene ) %>% gather(cell,counts,-Gene) %>% group_by(Gene,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) %>% select(colnames(counts[[n]])) %>% bind_rows(counts[[n]] %>% filter(!Gene %in% dup.gene ) ) %>% tibble::column_to_rownames("Gene"))[GI,]
  metas[[n]] <- metas[[n]] %>% left_join(data.frame(cell=colnames(counts[[n]]),mt.perc=colSums(counts[[n]][mt.gene,])/colSums(counts[[n]])) %>% tbl_df() ,by="cell")
}


counts.all <-do.call("bind_cols",counts) 
meta.all <- do.call("bind_rows",metas)

#' saving files
print("save output")
save(counts.all,meta.all,file=paste0("tmp_data/","/HB_coseq_smt2","/HB_coseq_smt2.all.counts.meta.Rdata"))
save(gtf.anno,mt.gene,ribo.gene,file=paste0("tmp_data/","/HB_coseq_smt2","/HB_coseq_smt2.gene.meta.Rdata"))


