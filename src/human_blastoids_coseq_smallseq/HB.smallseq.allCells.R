#' ---
#' title: UMAP for blastoids small-seq
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


#' human small meta
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)

sel.type <- c("miRNA")

meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.meta.filter.rds")) 
counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))
miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]


#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(EML=devTime,RNA_EML="None",small_EML="None") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.data.ob.umap.rds")) %>% select(cell,EML,RNA_EML,small_EML),by="cell")

#'miRNA only
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]


#'reference cells
ELC.cells <- meta.filter  %>% filter(EML %in% c("ELC")) %>% pull(cell)
TLC.cells <- meta.filter  %>% filter(EML %in% c("TLC")) %>% pull(cell)

#'
co.fm.list <- readRDS(paste0("tmp_data/",TD,"/HB.small.coseq_miRNA.split.fm.list.rds"))
temp.mk <- (co.fm.list$EML_mk %>% filter(p_val_adj < 0.05) %>% group_by(gene) %>% top_n(1,-1*p_val_adj) %>% split(.,.$cluster))[c("ELC","TLC")] %>% lapply(function(x){x$gene %>% head(9)})

if (file.exists(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
  sel.exp <- readRDS(paste0("tmp_data/",TD,"/HB.smallseq.miRNA.norm.rds"))
  miRNA.FM <- readRDS(paste0("tmp_data/",TD,"/","HB.smallseq.miRNA.miRNA.FM.mk.rds"))
}else{
  #' only small-seq ones
  temp.M <- meta.filter %>% filter(batch=="HB.smallseq") 
  temp.cells <- temp.M$cell
  temp.counts <- counts.miRNA.filter[,temp.M$cell]
  temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
  temp.sel.expG <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(temp.sel.expG)
  data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  temp.sce.sf <- data.frame(cell=colnames(temp.sce),sf=temp.sce$sizeFactor) %>% tbl_df()
  sel.exp <- scuttle::normalizeCounts(temp.sce)
  
  temp.norm <- sel.exp
  data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
  #ngene=length(temp.sel.expG);npc=25;
  #ngene=250;npc=25;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE,vars.to.regress = "embryo_batch") %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
  ngene=250;npc=25;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
  temp.M %>% group_by(devTime,batch) %>% summarise(nCell=n_distinct(cell)) %>% spread(batch,nCell)
  data.temp <- data.ob  %>% FindClusters( resolution = 0.4,verbose = FALSE) # fix
  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "embryo",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = (meta.filter %>% filter(Notes=="ICM-enriched") %>% pull(cell)))+NoLegend()+NoAxes()+ggtitle("ICM-enriched"),
    DimPlot(data.temp,cells.highlight = (meta.filter %>% filter(Notes=="TE-enriched") %>% pull(cell)))+NoLegend()+NoAxes()+ggtitle("TE-enriched")
    #DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
    #FeaturePlot(data.temp,"miRNA.nExpG",label=T)+NoLegend()+NoAxes()
  )
  
  temp.plot <- list()
  for ( n in unique(data.temp@meta.data$embryo)) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$embryo==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
  FeaturePlot(data.temp,temp.mk$ELC)
  FeaturePlot(data.temp,temp.mk$TLC)
  
  
  temp.exclude.cells <- meta.filter %>% filter(EML=="unknown") %>% pull(cell)
  #' all cells together
  temp.M <- meta.filter %>% filter(!cell %in% temp.exclude.cells)
  temp.cells <- temp.M$cell
  temp.counts <- counts.miRNA.filter[,temp.M$cell]
  temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
  data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  temp.sce.sf <- data.frame(cell=colnames(temp.sce),sf=temp.sce$sizeFactor) %>% tbl_df()
  sel.exp <- scuttle::normalizeCounts(temp.sce)
  temp.norm <- sel.exp
  data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
  #ngene=length(temp.sel.expG);npc=25;
  ngene=250;npc=25;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
  temp.M %>% group_by(devTime,batch) %>% summarise(nCell=n_distinct(cell)) %>% spread(batch,nCell)
  data.temp <- data.ob  %>% FindClusters( resolution = 1,verbose = FALSE) # fix
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = ELC.cells)+ggtitle("ELC")+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = TLC.cells)+ggtitle("TLC")+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    #DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  # clear batch effect
  
  #' run CCA integration
  temp.exclude.cells <- meta.filter %>% filter(EML=="unknown") %>% pull(cell)
  expG.set <- list()
  for (b in unique(meta.filter$batch  %>% unique() %>% as.vector())) {
    temp.cell <- meta.filter %>% filter(batch==b) %>% pull(cell) %>% setdiff(temp.exclude.cells)
    temp.counts <- counts.miRNA.filter[,temp.cell]
    expG.set[[b]] <- rownames(temp.counts)[rowSums(temp.counts[,temp.cell] >=1) >=2]
  }
  lapply(expG.set,length)
  sel.expG <- unlist(expG.set) %>% unique() %>% as.vector()
  
  sce.ob <- list()
  for (b in unique(meta.filter$batch  %>% unique() %>% as.vector())) {
    print(b)
    temp.M <- meta.filter %>% filter(batch==b) %>% filter(!cell %in% temp.exclude.cells)
    temp.counts <- counts.miRNA.filter[,temp.M$cell]
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  mBN.sce.ob <- multiBatchNorm(sce.ob$HB.smallseq,sce.ob$HB.coseq.small)
  names(mBN.sce.ob) <- c("HB.smallseq","HB.coseq.small")
  sel.exp <- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  
  temp.M <- meta.filter %>% filter(cell %in% colnames(sel.exp ))
  temp.counts <- counts.miRNA.filter[,temp.M$cell]
  temp.sel.expG <- rownames(sel.exp)  #%>% setdiff(chrX.miRNA.unique)
  data.merge <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.merge@assays$RNA$data <- as.matrix(sel.exp[temp.sel.expG,colnames(data.merge)])
  data.spt <- SplitObject(data.merge, split.by = "batch") %>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=200)})
  
  #k.anchor = 5 ; ;k.score = 30 ;
  set.seed(123)
  k.filter = 50 ;npc <- 25;ngene =250;k.weight = 50;sel.od <- c("HB.coseq.small","HB.smallseq")
  data.spt <- data.spt[sel.od]
  anchor.features <- SelectIntegrationFeatures(object.list = data.spt,nfeatures = ngene,verbose=F)
  #data.anchors <- FindIntegrationAnchors(data.spt,  anchor.features = anchor.features ,k.anchor = k.anchor, k.filter = k.filter, k.score = k.score,verbose=F)
  #data.integrated <- IntegrateData(anchorset = data.anchors,  k.weight = k.weight,verbose=F)
  data.anchors <- FindIntegrationAnchors(data.spt,  anchor.features = anchor.features,k.filter = k.filter,verbose=F)
  data.integrated <- IntegrateData(anchorset = data.anchors,  verbose=F,k.weight = k.weight )
  
  DefaultAssay(object = data.integrated) <- "integrated"
  data.ob  <- ScaleData(object = data.integrated, verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) 
  data.temp <- data.ob%>% FindClusters( resolution = 0.4,verbose = FALSE)
  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    #DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "small_EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = ELC.cells)+ggtitle("ELC")+NoLegend()+NoAxes(),
    DimPlot(data.temp,cells.highlight = TLC.cells)+ggtitle("TLC")+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    #DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,devTime,embryo,miRNA.nExpG,miRNA.UMI,batch,EML,small_EML,RNA_EML)) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")  
  data.ob.umap <- data.ob.umap %>% filter(batch=="HB.coseq.small") %>% mutate(EML=small_EML) %>% bind_rows(data.ob.umap %>% filter(batch=="HB.smallseq") %>% mutate(EML=recode(SC,"C0"="TLC","C1"="ELC")))
  data.ob.umap %>% group_by(devTime,EML,SC,batch) %>% summarise(nCell=n_distinct(cell))
  data.ob.umap %>% group_by(devTime,EML,batch) %>% summarise(nCell=n_distinct(cell))
  data.ob@meta.data$EML <-  (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
  data.ob@meta.data$SC <-  (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"SC"]
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))

}



