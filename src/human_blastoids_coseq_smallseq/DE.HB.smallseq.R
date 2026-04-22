#' ---
#' title: DE between HB ELC vs TLC (small-seq )
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
  library(batchelor)
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


raw.DEG.stat=matrix(nrow=2,ncol=1)
rownames(raw.DEG.stat) <- c("up_regulated","down_regulated")
colnames(raw.DEG.stat) <- c(
  paste("ELC","vs","TLC",sep="_")
)

savefile <- paste0("tmp_data/",TD,"/HB.smallseq.only.DEG.out.Rdata")
DEG.results.list <- list()
DEG.stat.list <- list()
sel.exp.list <- list()

sel.type <- c("miRNA","snoRNA","tRNA","snRNA","piRNA","rRNA")


if (file.exists(savefile)){
  load(savefile,verbose = T)
}else{
  #' human small meta
  load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)
  
  
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))
  
  table(duplicated(meta.filter$cell))
  table(meta.filter$batch)
  
  meta.filter <- meta.filter %>% filter(!batch %in% c("HB.coseq.small")) 
  cells.list <- meta.filter  %>% filter(devTime!="E5") %>% mutate(SID=EML) %>% split(.$SID) %>% lapply(function(x){x$cell})
  
  
  
  for ( st in sel.type) {
    DEG.stat <- raw.DEG.stat
    miRNA.ID <- trans.anno$mature %>% filter(type==st) %>% pull(ID) %>% intersect(counts.filter$mature$ID)
    chr.miRNA.ID <- miRNA.ID  %>% intersect(chr.smallRNA.id) 
    counts.miRNA.filter <- counts.filter$mature %>% filter(ID %in% chr.miRNA.ID)   %>% tibble::column_to_rownames("ID")
    counts.miRNA.filter <- counts.miRNA.filter[,meta.filter$cell]       
    
    rownames(counts.miRNA.filter) <- gsub("_","-",rownames(counts.miRNA.filter))
    
    expG.set <- list()
    for (b in unique(meta.filter$batch  %>% unique() %>% as.vector())) {
      temp.cell <- meta.filter %>% filter(batch==b) %>% pull(cell)
      temp.counts <- counts.miRNA.filter[,temp.cell]
      expG.set[[b]] <- rownames(temp.counts)[rowSums(temp.counts[,temp.cell] >=1) >=2]
    }
    lapply(expG.set,length)
    sel.expG <- unlist(expG.set) %>% unique() %>% as.vector()
    
    
    temp.M <- meta.filter 
    temp.counts <- counts.miRNA.filter[,temp.M$cell]
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sel.exp <- scuttle::normalizeCounts(temp.sce)
    
    sel.exp.list[[st]] <-  sel.exp 
    
    temp.M <- meta.filter %>% filter(cell %in% colnames(sel.exp ))
    temp.counts <- counts.miRNA.filter[,temp.M$cell]
    temp.sel.expG <- rownames(sel.exp)  #%>% setdiff(chrX.miRNA.unique)
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    data.ob@assays$RNA$data <- as.matrix(sel.exp[temp.sel.expG,colnames(data.ob)])
    
    
    DEG.results <- list()
    for (n in colnames(DEG.stat)) {
      DEG.results[[n]] <- n
    }
    
    for (n in names(DEG.results)) {
      print(n)
      DEG.results[[n]] <- FunDEG_custom(data.ob,cells.list,n)
    }
    
    for (n in names(DEG.results)) {
      print(n)
      DEG.stat["up_regulated",n] <- nrow(DEG.results[[n]][["DEG.result.up"]])
      DEG.stat["down_regulated",n]  <- nrow(DEG.results[[n]][["DEG.result.down"]])
      print(paste(DEG.stat["up_regulated",n],DEG.stat["down_regulated",n]))
    }
    DEG.results.list[[st]] <- DEG.results
    DEG.stat.list[[st]] <- DEG.stat
  }
  save(DEG.results.list,DEG.stat.list,sel.exp.list,file=savefile)
}




