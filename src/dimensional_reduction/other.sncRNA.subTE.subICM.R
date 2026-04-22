#' ---
#' title: (test version) sub ICM and sub TE cells based on other type ncRNA
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
  library(slingshot)
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
TD="Nov_2024"

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter


para_cal=TRUE
if (para_cal) {
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  numCores <- 11
  registerDoParallel(numCores)
}


#' loading gene meta data
load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)


sel.exp.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".norm.list.rds"))
sel.exp.list$miRNA <-  readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))
counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
# all cells  annotation
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds"))%>% mutate(sub_EML=EML) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(EML=devTime,RNA_EML="None",small_EML="None") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% select(cell,EML,RNA_EML,small_EML),by="cell") %>% mutate(sub_EML=EML)
#' update the full cell annotation
meta.filter <- meta.filter  %>% rows_update(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")


if (file.exists(paste0("tmp_data/",TD,"/","other.sncRNA",".sub.data.ob.list.rds"))) {
  data.ncRNA.sub.ob.list <- readRDS(paste0("tmp_data/",TD,"/","other.sncRNA",".sub.data.ob.list.rds"))
  data.ncRNA.sub.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","other.sncRNA",".sub.data.ob.umap.list.rds"))
  data.ncRNA.sub.mk.list <- readRDS(paste0("tmp_data/",TD,"/","other.sncRNA",".sub.ncRNA.sub.mk.list.rds"))
}else{
  data.ncRNA.sub.ob.list <- list()
  data.ncRNA.sub.ob.umap <- list()
  data.ncRNA.sub.mk.list <- list()
  for ( sel.type in c("snoRNA","snRNA","tRNA")) {
    sncRNA.ID <- trans.anno$mature %>% filter(type==sel.type) %>% pull(ID) %>% intersect(counts.filter$mature$ID)
    chr.sncRNA.ID <- trans.anno$mature %>% filter(type==sel.type) %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
    counts.sncRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.sncRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]
    sel.exp <-  sel.exp.list[[sel.type]]
    
    data.ob.list <- list()
    data.sub.mk.list <- list()
    #'  for ICM ones
    temp.M <- meta.filter %>% filter(EML=="ICM") %>% filter(batch %in% c("batch1","batch2","batch4"))
    temp.cells <- temp.M$cell
    temp.counts <- counts.sncRNA.filter[,temp.M$cell]
    temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    temp.norm <- sel.exp
    data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
    #ngene=length(temp.sel.expG);npc=25;
    ngene=100;npc=10;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
    temp.M %>% group_by(devTime,batch) %>% summarise(nCell=n_distinct(cell)) %>% spread(batch,nCell)
    
    data.temp <- data.ob  %>% FindClusters( resolution = 0.8,verbose = FALSE) # fix
    cowplot::plot_grid(
      DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "sub_EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
    )
    data.ob.list[["ICM"]] <- data.temp
    Idents(data.ob.list[["ICM"]]) <- factor( data.ob.list[["ICM"]]@meta.data$sub_EML)
    data.sub.mk.list[["ICM"]] <- FindAllMarkers(data.ob.list[["ICM"]],only.pos = T,logfc.threshold=0.25,min.pct=1/3) %>% tbl_df() %>% filter(p_val_adj < 0.05) %>% inner_join(log1p(AverageExpression(data.ob.list[["ICM"]])$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()  %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
    
    data.ob.umap <- data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()%>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")%>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% mutate(IT="ICM_sub",type=sel.type)
    
    #'  for TE ones
    temp.M <- meta.filter %>% filter(EML=="TE") %>% filter(batch %in% c("batch1","batch2","batch4"))
    temp.cells <- temp.M$cell
    temp.counts <- counts.sncRNA.filter[,temp.M$cell]
    temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    temp.norm <- sel.exp
    data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
    #ngene=length(temp.sel.expG);npc=25;
    ngene=100;npc=10;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
    
    data.temp <- data.ob  %>% FindClusters( resolution = 0.8,verbose = FALSE) # fix
    cowplot::plot_grid(
      DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "sub_EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
    )
    data.ob.list[["TE"]] <- data.temp
    Idents(data.ob.list[["TE"]]) <- factor( data.ob.list[["TE"]]@meta.data$sub_EML)
    data.sub.mk.list[["TE"]] <- FindAllMarkers(data.ob.list[["TE"]],only.pos = T,logfc.threshold=0.25,min.pct=1/3) %>% tbl_df() %>% filter(p_val_adj < 0.05) %>% inner_join(log1p(AverageExpression(data.ob.list[["TE"]])$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()  %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
    
    data.ob.umap <- data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()%>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")%>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% mutate(IT="TE_sub",type=sel.type) %>% bind_rows(data.ob.umap )
    
    data.ncRNA.sub.ob.list[[sel.type]] <- data.ob.list
    data.ncRNA.sub.ob.umap[[sel.type]] <- data.ob.umap
    data.ncRNA.sub.mk.list[[sel.type]]<- data.sub.mk.list
  }
  saveRDS(data.ncRNA.sub.ob.list,paste0("tmp_data/",TD,"/","other.sncRNA",".sub.data.ob.list.rds"))
  saveRDS(data.ncRNA.sub.ob.umap,paste0("tmp_data/",TD,"/","other.sncRNA",".sub.data.ob.umap.list.rds"))
  saveRDS(data.ncRNA.sub.mk.list,paste0("tmp_data/",TD,"/","other.sncRNA",".sub.ncRNA.sub.mk.list.rds"))
}
