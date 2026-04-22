#' ---
#' title: UMAP based on other ncRNA for all cells
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

# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=ifelse(EML=="prelineage" & devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(EML=="prelineage" & !devTime %in% c("2C","4C"),"L8CM",EML))  %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

#' update the full cell annotation
meta.filter <- meta.filter  %>% left_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell") ### need to be left_join

#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap,by="cell")
meta.filter <- meta.filter %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))

#' fix the "_" issues
trans.anno$mature <- trans.anno$mature %>% mutate(ID=gsub("_","-",ID))
trans.anno$prec <- trans.anno$prec %>% mutate(ID=gsub("_","-",ID))
counts.filter$mature <- counts.filter$mature %>% mutate(ID=gsub("_","-",ID))
chr.smallRNA.id <- gsub("_","-",chr.smallRNA.id)

if (file.exists(paste0("tmp_data/",TD,"/","allCells.other.sncRNA.ave.exp.list.rds"))) {
  sncRNA.ave.exp.list <- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA.ave.exp.list.rds"))
  data.ob.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.list.rds"))
  data.ob.umap.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.umap.list.rds"))
  sel.exp.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".norm.list.rds"))
  #sncRNA.FM.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA.FM.mk.list.rds"))
  #sncRNA.mk.out.list <- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA.mk.list.rds"))
}else{
  sncRNA.ave.exp.list <- list()
  data.ob.list <- list()
  data.ob.umap.list <- list()
  sel.exp.list <- list()
 
  for ( sel.type in c("snoRNA","snRNA","tRNA","rRNA","piRNA")) {
    sncRNA.ID <- trans.anno$mature %>% filter(type==sel.type) %>% pull(ID) %>% intersect(counts.filter$mature$ID)
    chr.sncRNA.ID <- trans.anno$mature %>% filter(type==sel.type) %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
    counts.sncRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.sncRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]
    
    
    
    #'ds size
    counts.sncRNA.filter <- (counts.filter$mature %>% filter(ID %in% sncRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]
    
    #'  IT
    temp.M <- meta.filter %>% filter(EML!="unknown") %>% filter(!batch %in% c("Split"))
    temp.cells <- temp.M$cell
    temp.counts <- counts.sncRNA.filter[,temp.M$cell]
    temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
    temp.sce.sf <- data.frame(cell=colnames(temp.sce),sf=temp.sce$sizeFactor) %>% tbl_df()
    sel.exp <- scuttle::normalizeCounts(temp.sce)
    temp.norm <- sel.exp
    data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
    #ngene=length(temp.sel.expG);npc=25;
    if (sel.type !="piRNA") {
      ngene=250;npc=25;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
      data.temp <- data.ob  %>% FindClusters( resolution = 0.8,verbose = FALSE) 
    }else{
      ngene=2000;npc=25;data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE)
      data.temp <- data.ob  %>% FindClusters( resolution = 0.8,verbose = FALSE) 
    }
    # fix
    cowplot::plot_grid(
      DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "sub_EML",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
      DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
    )
    
    data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,devTime,stage,paste0(sel.type,".nExpG"),paste0(sel.type,".UMI"),batch,EML,sub_EML,RNA_EML,small_EML)) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")  
    data.ob@meta.data$SC <-  (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"SC"]
    
    #' create average expression
    data.deg <- subset(data.ob,cell=(data.ob.umap %>% filter(batch!="Split" & EML!="unknown") %>% pull(cell)))
    Idents(data.deg) <- as.factor((data.deg@meta.data$EML))
    sncRNA.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() # exclude the split batch ones
    
    sncRNA.ave.exp.list[[sel.type]] <- sncRNA.ave.exp
    data.ob.list[[sel.type]] <- data.ob
    data.ob.umap.list[[sel.type]] <- data.ob.umap
    sel.exp.list[[sel.type]] <- sel.exp
    
  }
  saveRDS(sncRNA.ave.exp.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA.ave.exp.list.rds"))
  saveRDS(data.ob.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.list.rds"))
  saveRDS(data.ob.umap.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.umap.list.rds"))
  saveRDS(sel.exp.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".norm.list.rds"))
  #saveRDS(sncRNA.FM.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA.FM.mk.list.rds"))
  #saveRDS(sncRNA.mk.out.list,file=paste0("tmp_data/",TD,"/","allCells.other.sncRNA.mk.list.rds"))
  #saveRDS(sncRNA.em.mk.out.list,file=paste0("tmp_data/",TD,"/","only.embryonic.cells.other.sncRNA.mk.list.rds"))
}

