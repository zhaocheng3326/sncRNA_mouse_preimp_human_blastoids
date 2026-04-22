#' ---
#' title: "other nc UMAP for HB_coseq small"
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

#' human small meta
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)


#meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.meta.filter.rds")) 
meta.filter <-  readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.umap.rds")) %>% rename(SC=seurat_clusters) %>% select(cell,embryo,devTime,SC,EML,RNA_EML,small_EML)
counts.filter <- (readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))$mature%>% tibble::column_to_rownames("ID"))[,meta.filter$cell]

#' type of ncRNA
sel.type <- c("tRNA","snoRNA","piRNA","snRNA",'rRNA')
data.ob.list <- list()

rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))
trans.anno$mature$ID <- gsub("_","-",trans.anno$mature$ID)
data.ob.umap.list <- list()
data.ob.list <- list()
DE.ncRNA.list <- list()
if (file.exists(paste0("tmp_data/",TD,"/HB.coseq.small.otherDE.rds"))) {
  DE.ncRNA.list <- readRDS(paste0("tmp_data/",TD,"/HB.coseq.small.otherDE.rds"))
  data.ob.umap.list <- readRDS(paste0("tmp_data/",TD,"/HB.coseq.small.otherncRNA.umap.rds"))
}else{
  for (s in sel.type) { ## snoRNA also is a marker
    npc <- 10
    temp.M <- meta.filter #%>% filter(cell %in% colnames(sel.exp ))
    temp.sel.expG <- trans.anno$mature %>% filter(type==s) %>% pull(ID) %>% intersect(rownames(counts.filter))
    temp.counts <- counts.filter[temp.sel.expG,temp.M$cell]
    temp.M <- temp.M[colSums(temp.counts>0) > 30,]
    temp.cells <- temp.M$cell
    temp.sel.expG<-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
    temp.counts <- temp.counts[temp.sel.expG,temp.cells]
    
    data.temp <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    temp.norm <- scuttle::normalizeCounts(temp.sce)
    data.temp@assays$RNA$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.temp)])
    if (s=="piRNA") {
      data.temp <- data.temp %>% FindVariableFeatures(verbose=F,nfeatures=500) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
    }else{
      #data.temp <- data.temp %>% FindVariableFeatures(verbose=F,nfeatures=length(temp.sel.expG)) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1>
      data.temp <- data.temp %>% FindVariableFeatures(verbose=F,nfeatures=150) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
    }
    
    DimPlot(data.temp,group.by = "EML")+ggtitle(s)
    data.ob.list[[s]] <- data.temp
    data.deg.temp <- subset(data.temp,cell=(temp.M %>% filter(EML!="other") %>% pull(cell)))
    Idents(data.deg.temp) <- factor(data.deg.temp@meta.data$EML)
    DE.ncRNA.list[[s]] <- FindAllMarkers(data.deg.temp,verbose = F,only.pos = T)  %>% tbl_df()

    data.ob.umap.list[[s]] <- Embeddings(data.ob.list[[s]],reduction='umap') %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()
  }
  #' saving object
  saveRDS(DE.ncRNA.list ,paste0("tmp_data/",TD,"/HB.coseq.small.otherDE.rds"))
  saveRDS(data.ob.umap.list ,paste0("tmp_data/",TD,"/HB.coseq.small.otherncRNA.umap.rds"))
}
