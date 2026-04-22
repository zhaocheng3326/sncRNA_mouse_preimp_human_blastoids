#' ---
#' title: UMAP based on miRNA integration
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



if (file.exists(paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.data.ob.umap.rds"))
  sel.exp <- readRDS(paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.norm.rds"))
}else{
  HB.meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
  HB.counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))
  
  
  HE.meta.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.coord.check.rds")%>% mutate(rename_EML=ifelse(EML %in% c("EarlyTE","mural","polar"),"TE",EML)) %>% mutate(EML=rename_EML) %>% select(cell,embryo,devTime,batch,EML)
  HE.counts.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")
  
  
  meta.filter <- HB.meta.filter %>% select(cell,embryo,batch,EML,devTime) %>% bind_rows(HE.meta.filter %>% select(cell,embryo,batch,EML,devTime) ) %>% filter(!batch %in% c("HB.coseq.small","HumanSplit")) 
  table(duplicated(meta.filter$cell))
  table(meta.filter$batch)
  
  miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(HB.counts.filter$mature$ID)%>% intersect(HE.counts.filter$mature$ID)
  chr.miRNA.ID <- miRNA.ID  %>% intersect(chr.smallRNA.id) 
  counts.miRNA.filter <- HB.counts.filter$mature %>% filter(ID %in% chr.miRNA.ID)  %>% inner_join(HE.counts.filter$mature %>% filter(ID %in% chr.miRNA.ID),by="ID")  %>% tibble::column_to_rownames("ID")
  counts.miRNA.filter <- counts.miRNA.filter[,meta.filter$cell]                      
  
  
  expG.set <- list()
  for (b in unique(meta.filter$batch  %>% unique() %>% as.vector())) {
    temp.cell <- meta.filter %>% filter(batch==b) %>% pull(cell) 
    temp.counts <- counts.miRNA.filter[,temp.cell]
    expG.set[[b]] <- rownames(temp.counts)[rowSums(temp.counts[,temp.cell] >=1) >=2]
  }
  lapply(expG.set,length)
  sel.expG <- unlist(expG.set) %>% unique() %>% as.vector()
  
  sce.ob <- list()
  for (b in unique(meta.filter$batch  %>% unique() %>% as.vector())) {
    print(b)
    temp.M <- meta.filter %>% filter(batch==b) 
    temp.counts <- counts.miRNA.filter[,temp.M$cell]
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  mBN.sce.ob <- multiBatchNorm(sce.ob$SB,sce.ob$B1,sce.ob$TB,sce.ob$A1,sce.ob$HB.smallseq)
  names(mBN.sce.ob) <- c("SB","B1","TB","A1","HB.smallseq")
  sel.exp <- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  
  
  temp.exclude.cells <- meta.filter %>% filter(EML=="unclassified" ) %>% pull(cell)
  temp.M <- meta.filter %>% filter(cell %in% colnames(sel.exp )) %>% filter(!cell %in% temp.exclude.cells)
  temp.counts <- counts.miRNA.filter[,temp.M$cell]
  temp.sel.expG <- rownames(sel.exp)  
  data.merge <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.merge@assays$RNA$data <- as.matrix(sel.exp[temp.sel.expG,colnames(data.merge)])
  data.spt <- SplitObject(data.merge, split.by = "batch") %>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=200)})
  
  #k.anchor = 5 ; ;k.score = 30 ;
  set.seed(123)
  k.filter = 50 ;npc <- 20;ngene =250;k.weight = 50;sel.od <- c("SB","B1","TB","A1","HB.smallseq")
  data.spt <- data.spt[sel.od]
  anchor.features <- SelectIntegrationFeatures(object.list = data.spt,nfeatures = ngene,verbose=F)
  #data.anchors <- FindIntegrationAnchors(data.spt,  anchor.features = anchor.features ,k.anchor = k.anchor, k.filter = k.filter, k.score = k.score,verbose=F)
  #data.integrated <- IntegrateData(anchorset = data.anchors,  k.weight = k.weight,verbose=F)
  data.anchors <- FindIntegrationAnchors(data.spt,  anchor.features = anchor.features,k.filter = k.filter,verbose=F)
  data.integrated <- IntegrateData(anchorset = data.anchors,  verbose=F,k.weight = k.weight )
  
  DefaultAssay(object = data.integrated) <- "integrated"
  data.ob  <- ScaleData(object = data.integrated, verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) 
  data.temp <- data.ob%>% FindClusters( resolution = 0.4,verbose = FALSE)
  
  data.ob.umap <- temp.M %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell")
  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
    #DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  temp.plot <- list()
  for ( n in unique(HB.meta.filter %>% filter(batch=="HB.smallseq") %>% pull(embryo))) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[ data.temp@meta.data$embryo==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.data.ob.umap.rds"))
  saveRDS(sel.exp,paste0("tmp_data/",TD,"/HB.HE.CCA.miRNA.norm.rds"))
  
  
  #' get log2FC
  # embryo
  #' generate the average expression
  data.human.umap <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.coord.check.rds")
  human.norm.exp <- readRDS(paste0("~/My_project/Extra_smncRNA/tmp_data/","Apr_2023","/","miRNA",".norm.rds"))$mBN #%>% tibble::column_to_rownames("ID")

  
  data.temp.umap <- data.human.umap %>% filter(batch!="HumanSplit") %>% mutate(SID=recode(EML,"EarlyTE"="TE","mural"="TE","polar"="TE","E3"="L8CM","E4"="L8CM")) %>% filter(SID %in% c("ICM","TE"))
  temp.counts <- (readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")$mature %>% tibble::column_to_rownames("ID"))[rownames(human.norm.exp),data.temp.umap$cell]
  
  data.sub.temp <- CreateSeuratObject(temp.counts[rownames(human.norm.exp),], meta.data = (data.temp.umap %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.sub.temp@assays$RNA$data <- as.sparse(human.norm.exp[,colnames(data.sub.temp@assays$RNA$data)])
  
  data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  Idents(data.sub.temp) <- factor(data.sub.temp@meta.data$SID)
  DefaultAssay(data.sub.temp) <- "RNA"
  temp.ICMvsTE.log2FC <- FindMarkers(data.sub.temp,ident.1="ICM",ident.2="TE",logfc.threshold=0,min.pct =0) %>% tibble::rownames_to_column("gene") %>% tbl_df()
  
  #' blastoids
  data.human.umap <-  readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
  data.temp <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.rds"))
 
  data.temp.umap <- data.human.umap %>% filter(batch!="HB.coseq.small") %>% mutate(SID=EML) %>% filter(SID %in% c("ELC","TLC"))
  data.sub.temp <- subset(data.temp,cells=data.temp.umap$cell)
  data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  Idents(data.sub.temp) <- factor(data.sub.temp@meta.data$SID)
  DefaultAssay(data.sub.temp) <- "RNA"
  temp.ELCvsTLC.log2FC <- FindMarkers(data.sub.temp,ident.1="ELC",ident.2="TLC",logfc.threshold=0,min.pct =0) %>% tibble::rownames_to_column("gene") %>% tbl_df()
  
  IvM.EvT.fm <- temp.ICMvsTE.log2FC %>%mutate(type="ICMvsTE") %>% bind_rows(temp.ELCvsTLC.log2FC%>%mutate(type="ELCvsTLC") )
  saveRDS(IvM.EvT.fm,paste0("tmp_data/",TD,"/HB.HE.IvM.EvT.fm.rds"))
  
   
  #' Average expression (not.used)
  #data.temp.umap <- data.ob.umap %>% filter(EML %in% c("ICM","TE","ELC","TLC")) %>% mutate(SID=EML)
  #data.sub.temp <- subset(data.ob,cells=data.temp.umap$cell)
  #DefaultAssay(data.sub.temp) <- "RNA"
  #data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  #Idents(data.sub.temp) <- factor(data.sub.temp@meta.data$SID)
  #data.HE.HB.ave.exp <- AverageExpression(data.sub.temp,group.by = "SID")$RNA %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% gather(cluster,ave_exp,-gene) %>% tbl_df() %>% mutate(species="human") %>% mutate(cluster=gsub("-","_",cluster)) %>% rename(SID=cluster)
  #saveRDS(data.HE.HB.ave.exp,paste0("tmp_data/",TD,"/HB.HE.miRNA.ave.exp.rds"))
  
}


data.temp <- data.ob%>% FindClusters( resolution = 0.4,verbose = FALSE)
cowplot::plot_grid(
  DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "batch",label=T)+NoLegend()+NoAxes(),
  #DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
  DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
)
temp.plot <- list()
for ( n in unique(data.ob.umap  %>% filter(batch=="HB.smallseq") %>% pull(embryo))) {
  temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[ data.temp@meta.data$embryo==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
}
print(cowplot::plot_grid(plotlist=temp.plot))



  