#' ---
#' title: "QC for HB_coseq small"
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

sel.type <- c("miRNA")

meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.meta.filter.rds")) 
counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))

miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]


#' running
sel.type <- c("miRNA")
zs.limit <- 2
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)


#'loading data (for HB coseq small only)
meta.filter <- meta.filter %>% filter(batch=="HB.coseq.small")
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]

#' loading coseq smt2 cell annotation
HB.Split.meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smt2.updated.data.ob.umap.rds"))

#'reference cells
ELC.cells <- HB.Split.meta.filter  %>% filter(EML %in% c("ELC")) %>% pull(cell)
TLC.cells <- HB.Split.meta.filter  %>% filter(EML %in% c("TLC")) %>% pull(cell)


temp.M <- meta.filter %>% filter(batch=="HB.coseq.small") 
temp.counts <- counts.miRNA.filter[,temp.M$cell]
temp.M <- temp.M %>% left_join(HB.Split.meta.filter %>% select(cell,EML),by="cell")  %>% mutate(EML=ifelse(is.na(EML),"None",EML))
temp.cells <- temp.M$cell
temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]


if (file.exists(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.umap.rds"))
  mk.list <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.split.mk.list.rds"))
  fm.list <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.split.fm.list.rds"))
  sel.exp <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.norm.rds"))
}else{
  # npc <- 25;ngene <- length(temp.sel.expG)
  npc <- 10;ngene <- 150
  data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
  temp.norm <- scuttle::normalizeCounts(temp.sce)
  data.ob@assays$RNA@layers$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.ob)])
  data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
  # data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
  
  data.temp <- data.ob  %>% FindClusters( resolution = 0.6,verbose = FALSE)
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = ELC.cells)+ggtitle("ELC")+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = TLC.cells)+ggtitle("TLC")+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% mutate(RNA_EML=EML)%>% select(c(cell,embryo,devTime,RNA_EML)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp)))) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_15),by="cell")
  
  data.ob.umap <- data.ob.umap %>% mutate(small_EML=recode(seurat_clusters,"C0"="TLC","C3"="TLC","C1"="ELC","C2"="ELC"))  %>% mutate(EML=ifelse(RNA_EML==small_EML | RNA_EML=="None" ,small_EML,"unknown"))
 
  data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
  data.ob@meta.data$small_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"small_EML"]
  data.ob@meta.data$RNA_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"RNA_EML"]
  data.ob@meta.data$seurat_clusters <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"seurat_clusters"]
  
  
  meta.filter.updataed <- meta.filter %>% left_join(data.ob.umap %>% select(cell,EML,RNA_EML,small_EML),by="cell")
  #' save object
  saveRDS(data.ob,paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.umap.rds"))
  saveRDS(meta.filter.updataed,paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.updated.meta.rds"))
  saveRDS(temp.norm,paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.norm.rds"))
  sel.exp <- temp.norm
  
}


