#' ---
#' title: "Define the lineages of split-seq cells" 
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


MSplit.counts.filter<- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.counts.filter.rds"))
MSplit.meta.filter <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.meta.filter.rds")) %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))%>% mutate(pj=recode(pj,"co_smt2_extra"="co_smt2"))


#' co-smt2 dataset only
temp.M <- MSplit.meta.filter %>% filter(pj=="co_smt2")
temp.sel.expG <- rownames(MSplit.counts.filter)[rowSums(MSplit.counts.filter[,c(temp.M$cell)] >=1) >2]

nGene=2000;npc=25; data.temp <- CreateSeuratObject(MSplit.counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:npc,verbose=F) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) %>%  FindClusters(resolution = 0.8,verbose = FALSE)  #,vars.to.regress=c("nGene") 

cowplot::plot_grid(
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
  #DimPlot(data.temp,group.by=c("EML"),label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by=c("pj"),label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by=c("stage"),label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by=c("devTime"),label=T)+NoAxes()+NoLegend()
  #FeaturePlot(data.temp,c("CCR7"))+NoAxes()+NoLegend()
)
#' save coseq smt2 expression
coseq.gene.exp <- data.temp@assays$RNA$data %>% as.data.frame()
saveRDS(coseq.gene.exp,paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))

#' loading Posfai et al., 2017
meta_ref <- readRDS(paste0("~/My_project/Gpig_scRNA/tmp_data/","Oct_2022","/mouse.meta.filter.rds")) %>% filter(pj=="Posfai_2017")
remove.cells <- readRDS(paste0("~/My_project/Gpig_scRNA/tmp_data/","Oct_2022","/Posfai_2017.remove.cells.rds"))
meta_ref <- meta_ref %>% filter(!cell %in% remove.cells)

count_ref <- readRDS(paste0("~/My_project/Gpig_scRNA/tmp_data/","Oct_2022","/mouse.counts.filter.rds"))[,meta_ref$cell] 

table(rownames(count_ref) %in% rownames(MSplit.counts.filter))
table(rownames(count_ref) == rownames(MSplit.counts.filter))

GI <- rownames(count_ref) %>% intersect(rownames(MSplit.counts.filter))

meta.filter <- MSplit.meta.filter %>% bind_rows(meta_ref) 
counts.filter <- cbind(MSplit.counts.filter[GI,],count_ref[GI,])
sel.mk <- c("Pdgfra","Sox17","Bmp2","Gata4","Foxa2","Gata6","Dppa1","Cdx2","EOMES","Krt8","Krt18","Nanog","Sox2","Upp1","Spp1")

#' loading lineage marker genes from published paper
if (file.exists(paste0("tmp_data/",TD,"/mouse_coseq_smt2.data.it.ob.rds"))) {
  data.cca.temp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.data.it.ob.rds"))
  MSplit.meta.filter <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds"))
  data.cca.umap <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.data.ob.umap.rds"))
}else{
  
  #' loading counts
  
  
  #' prepare for integration
  temp.M <- meta.filter
  temp.sel.cell <- temp.M$cell
  sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >2] 
  data.merge <- CreateSeuratObject(counts.filter[sel.expG,temp.sel.cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  data.spt <-  data.spt[c("Posfai_2017","co_smt2")]
  
  
  #' cca 
  #k.weight = 50;nGene=1500;pc=25;k.filter=50 #k.anchor , k.score are the default values#k.anchor = 5;k.score = 30;
  nGene=2000;nPC=25;k.filter=100
  sel.features <- SelectIntegrationFeatures(object.list = data.spt, nfeatures = nGene,verbose=F)
  data.merge.anchors <- FindIntegrationAnchors(object.list = data.spt,  anchor.features = sel.features, verbose = FALSE, k.filter = k.filter)
  data.cca.temp <- IntegrateData(anchorset = data.merge.anchors, verbose = FALSE)%>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE,npcs = 50)%>%  RunUMAP(dims = 1:nPC,verbose = FALSE)%>% FindNeighbors( dims = 1:nPC,verbose = FALSE)
  data.temp <- data.cca.temp%>%  FindClusters(resolution = 0.4,verbose = FALSE)
  #data.temp <- data.cca.temp%>%  FindClusters(resolution = 0.4,verbose = FALSE,algorithm = 2)
  data.sub.temp <- subset(data.temp,cell=(meta.filter %>% filter(pj=="co_smt2") %>% pull(cell)))
  DefaultAssay(data.sub.temp) <- "RNA"
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("EML"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("pj"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("stage"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.sub.temp,group.by=c("devTime"),label=T)+NoAxes()+NoLegend()
    #FeaturePlot(data.temp,c("CCR7"))+NoAxes()+NoLegend()
  )
  
  
  temp.M.updated <- MSplit.meta.filter %>% select(cell,devTime,pj) %>% filter(pj=="co_smt2")
  temp.M.updated$SC <- paste("C",Idents(data.sub.temp)[temp.M.updated$cell],sep="")
  temp.M.updated <- temp.M.updated %>% mutate(EML=recode(SC,"C0"="TE","C1"="ICM","C2"="prelineage","C3"="prelineage","C4"="prelineage") ) %>% mutate(EML=ifelse(devTime %in% c("64C") & EML=="prelineage","unknown",EML))%>% mutate(EML=ifelse(devTime %in% c("8C","2C","4C") & EML %in% c("ICM","TE"),"unknown",EML))
  temp.M.updated$EML %>% table()
  
  data.sub.temp <- subset(data.merge,cell=temp.M.updated$cell)
  data.sub.temp@meta.data$EML <- (temp.M.updated %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"EML"]
  VlnPlot(data.sub.temp,group.by = "EML",intersect(sel.mk,temp.sel.expG),ncol=7)
  #meta.filter <- meta.filter %>% rows_update(HSplit.meta.filter %>% select(cell,EML,big_EML),by="cell")
  
  MSplit.meta.filter <- MSplit.meta.filter %>% mutate(EML=paste0("L",devTime)) %>%  rows_update(temp.M.updated %>% select(cell,EML),by="cell")
  
  #DefaultAssay(data.cca.temp)="RNA"
  #saveRDS(data.mnn.temp,file=paste0("tmp_data/",TD,"/mouse_coseq_smt2.data.it.ob.rds"))
  saveRDS(data.cca.temp,file=paste0("tmp_data/",TD,"/mouse_coseq_smt2.data.it.ob.rds"))
  saveRDS(MSplit.meta.filter,paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds"))
  
  data.cca.umap <- Embeddings(data.cca.temp,reduction="umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% inner_join( data.cca.temp@meta.data %>% tibble::rownames_to_column("cell") %>% select(cell,devTime,pj,EML),by="cell") %>% tbl_df() %>% rows_update(MSplit.meta.filter %>% select(cell,EML),by="cell")
  #data.mnn.umap <- Embeddings(data.mnn.temp,reduction="umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% inner_join( data.mnn.temp@meta.data %>% tibble::rownames_to_column("cell") %>% select(cell,devTime,pj,EML),by="cell") %>% tbl_df() %>% rows_update(MSplit.meta.filter %>% select(cell,EML),by="cell")
  #saveRDS(data.mnn.umap,paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.data.ob.umap.rds"))
  saveRDS(data.cca.umap,paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.data.ob.umap.rds"))
  
  #" get the ICM vs TE DEGs
  #' co-smt2 dataset only
  temp.M <- data.cca.umap %>% filter(pj=="co_smt2") 
  temp.sel.expG <- rownames(MSplit.counts.filter)[rowSums(MSplit.counts.filter[,c(temp.M$cell)] >=1) >2]
  
  nGene=2000;npc=25; data.temp <- CreateSeuratObject(MSplit.counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:npc,verbose=F) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) %>%  FindClusters(resolution = 0.8,verbose = FALSE)  #,vars.to.regress=c("nGene") 
  
  Idents(data.temp) <- factor(data.temp@meta.data$EML)
  temp.DEG <- FindMarkers(data.temp,ident.1="ICM",ident.2="TE") %>% tibble::rownames_to_column("gene")%>% tbl_df()
  saveRDS(temp.DEG,paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.ICMvsTE.DEG.rds"))
  
  
}
