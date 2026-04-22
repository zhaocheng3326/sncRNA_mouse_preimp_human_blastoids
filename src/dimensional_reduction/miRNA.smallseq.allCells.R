#' ---
#' title: UMAP based on miRNA (coseq+smallseq)
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

meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)  %>% filter(batch!="batch5")
meta.filter <- meta.filter %>% filter(batch !="batch4") %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]

#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(EML=devTime,RNA_EML="None",small_EML="None") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% select(cell,EML,RNA_EML,small_EML),by="cell") %>% filter(EML!="unknown") # exclude incons cells from co-seq
meta.filter <- meta.filter %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))


#'miRNA only
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]


# loading norm exp 

if (file.exists(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))) {
   data.ob <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))
   data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds"))
   sel.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))
   miRNA.ave.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.ave.exp.rds"))
   #miRNA.FM <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.FM.mk.rds"))
   #miRNA.mk.out <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.mk.rds"))
}else{
  #' 
  temp.M <- meta.filter 
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
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,devTime,stage,miRNA.nExpG,miRNA.UMI,batch)) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% mutate(EML=ifelse(devTime %in% c("oocyte","sperm"),devTime,"none"))  %>% mutate(EML=ifelse(devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(devTime %in% c("8C") & SC %in% c("C2","C3","C4")   ,"L8CM",EML))  %>% mutate(EML=ifelse(EML=="none" & SC %in% c("C2"),"L8CM",EML)) %>% mutate(EML=ifelse(EML=="none" & SC %in% c("C0","C4"),"TE",EML)) %>% mutate(EML=ifelse(EML=="none" & SC=="C1","ICM",EML))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")  
  data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(devTime %in% c("8C") & EML %in% c("ICM","TE"),"unknown",EML)) %>% mutate(EML=ifelse(devTime %in% c("64C") & EML %in% c("L8CM"),"unknown",EML))
  data.ob.umap %>% group_by(devTime,EML,SC) %>% summarise(nCell=n_distinct(cell))
  data.ob@meta.data$EML <-  (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
  data.ob@meta.data$SC <-  (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"SC"]
  
  saveRDS(data.ob,file=paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))
  saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds"))
  saveRDS(sel.exp,file=paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))
  saveRDS(temp.sce.sf,paste0("tmp_data/",TD,"/","allCells.miRNA",".sf.rds"))
  
  #' create average expression
  data.deg <- subset(data.ob,cell=(data.ob.umap %>% filter(batch!="Split" & EML!="unknown") %>% pull(cell)))
  Idents(data.deg) <- as.factor((data.deg@meta.data$EML))
  miRNA.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() # exclude the split batch ones
  saveRDS(miRNA.ave.exp,file=paste0("tmp_data/",TD,"/","allCells.miRNA.ave.exp.rds"))
  
}

#' check the typical markers
FeaturePlot(data.ob,tp.mk )
data.temp <- data.ob  %>% FindClusters( resolution = 0.8,verbose = FALSE) # fix
cowplot::plot_grid(
  DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
  DimPlot(data.ob  %>% FindClusters( resolution = 1.1,verbose = FALSE),label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
)


print(
  data.ob.umap %>% group_by(stage,EML) %>% summarise(nCell=n_distinct(cell)) %>% ggplot()+geom_bar(mapping=aes(fill=stage,x=EML,y=nCell),stat="identity",position="dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))
)

# some extra check
data.temp <- data.ob  %>% FindClusters( resolution = 0.6,verbose = FALSE) # fix
cowplot::plot_grid(
  DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,cells.highlight = ICM.cells)+ggtitle("ICM")+NoLegend()+NoAxes(),
  DimPlot(data.temp,cells.highlight = TE.cells)+ggtitle("TE")+NoLegend()+NoAxes(),
  #DimPlot(data.temp,cells.highlight = PE.cells)+ggtitle("sus_PE")+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
)


