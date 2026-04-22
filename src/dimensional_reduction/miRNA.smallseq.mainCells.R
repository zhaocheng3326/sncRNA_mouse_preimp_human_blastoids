#' ---
#' title: UMAP based on miRNA used for traj inferring
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

sel.type <- c("miRNA")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch) 
meta.filter <- meta.filter %>% filter(batch !="batch4") %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]


#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(EML=devTime,RNA_EML="None",small_EML="None") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% select(cell,EML,RNA_EML,small_EML),by="cell")


#' update the full cell annotation 
meta.filter <- meta.filter  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% select(cell,EML),by="cell")


#'miRNA only
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]


if (file.exists(paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.umap.rds"))
}else{
  #' 
  temp.M <- meta.filter %>% filter(devTime!="sperm" & devTime !="oocyte") %>% filter(batch %in% c("batch1","batch2","batch4"))
  temp.cells <- temp.M$cell
  temp.counts <- counts.miRNA.filter[,temp.M$cell]
  temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
  data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  sel.exp <- scuttle::normalizeCounts(temp.sce)
  temp.norm <- sel.exp
  data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
  #ngene=length(temp.sel.expG);npc=25;
  ngene=100;npc=10;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
  temp.M %>% group_by(devTime,batch) %>% summarise(nCell=n_distinct(cell)) %>% spread(batch,nCell)
  data.ob@reductions$umap@cell.embeddings[,1] <- -1*data.ob@reductions$umap@cell.embeddings[,1] ## reverse the UMAP_1 coord
  data.temp <- data.ob  %>% FindClusters( resolution = 0.6,verbose = FALSE) # fix
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()%>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")
    
  saveRDS(data.ob,file=paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.rds"))
  saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.umap.rds"))
  
  
}
