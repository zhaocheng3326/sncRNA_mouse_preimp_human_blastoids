#' ---
#' title: check the small RNA part for split seq (based on other type of ncRNA)
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


load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)

meta.filter <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(SC=seurat_clusters) %>% select(cell,embryo,sex,devTime,stage,RNA_EML,small_EML,EML,SC) 
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

counts.filter <- (readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))$mature%>% tibble::column_to_rownames("ID"))[,meta.filter$cell]


sel.type <- c("tRNA","snoRNA","piRNA","snRNA",'rRNA') 
#' check sub length tRNA
#' check the RNA expression with specific lengthtRNA.WP.out <- readRDS(paste0("tmp_data/",td,"/tRNA.WP.out.rds"))
#tRNA.WP.out <- readRDS(paste0("tmp_data/","merge_batch_soft_link","/tRNA.WP.out.rds"))
#tRNA.exp <- list()
#tRNA.exp$l37 <- tRNA.WP.out %>% lapply(function(x) {x$tRNAL37.exp }) %>% do.call("bind_rows",.)#%>% spread(cell,umi) %>% replace(.,is.na(.),0)
#tRNA.exp$l74 <- tRNA.WP.out %>% lapply(function(x) {x$tRNAL74.exp }) %>% do.call("bind_rows",.)#%>% spread(cell,umi) %>% replace(.,is.na(.),0)
#tRNA.exp$l75 <- tRNA.WP.out %>% lapply(function(x) {x$tRNAL75.exp }) %>% do.call("bind_rows",.)#%>% spread(cell,umi) %>% replace(.,is.na(.),0)

data.ob.list <- list()

rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))
trans.anno$mature$ID <- gsub("_","-",trans.anno$mature$ID)
data.ob.umap.list <- list()
data.ob.list <- list()
DE.ncRNA.list <- list()
if (file.exists(paste0("tmp_data/",TD,"/Msmall.coseq.small.otherDE.rds"))) {
  DE.ncRNA.list <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq.small.otherDE.rds"))
  data.ob.list <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq.small.otherncRNA.ob.rds"))
  data.ob.umap.list <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq.small.otherncRNA.umap.rds"))
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
      #data.temp <- data.temp %>% FindVariableFeatures(verbose=F,nfeatures=length(temp.sel.expG)) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
      data.temp <- data.temp %>% FindVariableFeatures(verbose=F,nfeatures=150) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
    }
    
    DimPlot(data.temp,group.by = "EML")+ggtitle(s)
    data.ob.list[[s]] <- data.temp
    data.deg.temp <- subset(data.temp,cell=(temp.M %>% filter(EML!="other") %>% pull(cell)))
    Idents(data.deg.temp) <- factor(data.deg.temp@meta.data$EML)
    DE.ncRNA.list[[s]] <- FindAllMarkers(data.deg.temp,verbose = F,only.pos = T)  %>% tbl_df() 
    
    data.deg.temp <- subset(data.temp,cell=(temp.M %>% filter(EML!="other") %>% pull(cell)))
    Idents(data.deg.temp) <- factor(data.deg.temp@meta.data$EML)
    data.deg.temp <- RenameIdents(data.deg.temp,"ICM"="ICM_TE","TE"="ICM_TE") 
    DE.ncRNA.list[[s]] <- FindAllMarkers(data.deg.temp,verbose = F,only.pos = T)  %>% tbl_df() %>% filter(cluster=="ICM_TE") %>% bind_rows(DE.ncRNA.list[[s]])
    data.ob.umap.list[[s]] <- Embeddings(data.ob.list[[s]],reduction='umap') %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()
  }
  
  #' saving object
  saveRDS(DE.ncRNA.list ,paste0("tmp_data/",TD,"/Msmall.coseq.small.otherDE.rds"))
  saveRDS(data.ob.list ,paste0("tmp_data/",TD,"/Msmall.coseq.small.otherncRNA.ob.rds"))
  saveRDS(data.ob.umap.list ,paste0("tmp_data/",TD,"/Msmall.coseq.small.otherncRNA.umap.rds"))
}
