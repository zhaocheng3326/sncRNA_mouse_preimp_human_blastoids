#' ---
#' title: "HB UMAP"
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
  condaENV <- "/home/chenzh/miniconda3/envs/R4.0" 
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

counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smt2.counts.filter.rds"))
meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smt2.meta.filter.rds"))


human.mk.list.short <- list(EPI=c("NANOG","POU5F1","HMGN3"),PE=c("SOX17","GATA4","PDGFRA"),TE=c("GATA2","GATA3","GCM1"),Amnion=c("ISL1","GABRP"),ExEMes=c("LUM","FOXF1"),STB=c("CGB1","NCF4","TCAF2C","CGB2"))

#' reference tool prediction results
predict_out <- readRDS(paste0("tmp_data/",TD,"/HB_co_smt2_ref.pred.rds"))
predict_out$umap %>% APPJ_devTime_UMAP("HB_co_smt2")
predict_out$full.anno %>% rename(cell=query_cell)  %>% inner_join(meta.filter %>% select(cell,devTime),by="cell") %>% group_by(devTime,sub_pred_EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(sub_pred_EML,nCell)

data.ref.pred <- list()
data.ref.pred$anno <- predict_out$full.anno %>% rename(cell=query_cell) %>% full_join(meta.filter,by="cell")
data.ref.pred$umap <-  predict_out$umap


temp.M <- meta.filter %>% left_join(data.ref.pred$anno %>% select(cell,pred_EML,sub_pred_EML),by="cell")
temp.sel.expG <-  rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >2]

data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)

data.ob <- data.ob  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE,k.param = 15)  #"mt.perc",#%>% CellCycleScoring(s.features = s.genes, g2m.features =g2m.genes)  #vars.to.regress=c("mt.perc","S.Score", "G2M.Score")

#data.ob <- data.ob  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F,vars.to.regress=c("nGene"))%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE)  #"mt.perc",#%>% CellCycleScoring(s.features = s.genes, g2m.features =g2m.genes)  #vars.to.regress=c("mt.perc","S.Score", "G2M.Score")

data.temp <- data.ob %>% FindClusters(reso=0.8,verbose=F)
plot_grid(
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
  #DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by="devTime",label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
  DimPlot(data.temp,group.by="batch",label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"POU5F1")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"NANOG")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"SOX2")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"GATA2")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"GATA3")+NoAxes()+NoLegend(),
  #FeaturePlot(data.temp,"CDX2")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
  #FeaturePlot(data.temp,"ribo.perc")+NoAxes()+NoLegend(),
 # DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend(),
  DimPlot(data.temp,label=T,group.by = "pred_EML")+NoAxes()+NoLegend(),
  DimPlot(data.temp,label=T,group.by = "sub_pred_EML")+NoAxes()+NoLegend()
)


temp.plot <- list()
for ( n in unique(data.temp@meta.data$embryo)) {
  temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[ data.temp@meta.data$embryo==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
}
print(cowplot::plot_grid(plotlist=temp.plot))

c("GCM1","CYP19A1","CCR7","OVOL1","PAPOLA","GNA12")%>% FunFP_plot(data.temp,.) %>% cowplot::plot_grid(plotlist = .)
c("PDGFRA","SOX17","GATA4","BMP2") %>% FunFP_plot(data.temp,.) %>% cowplot::plot_grid(plotlist = .)
unlist(human.mk.list.short) %>% FunFP_plot(data.temp,.) %>% cowplot::plot_grid(plotlist = .)


counts.filter[c("CYP19A1","CCR7","OVOL1","GCM1"),meta.filter$cell] %>% tibble::rownames_to_column("gene") %>% gather(cell,ct,-gene) %>% tbl_df() %>% filter(ct > 0) %>% inner_join(meta.filter %>% select(cell,embryo),by="cell") %>% group_by(gene,embryo) %>% summarise(nCell=n_distinct(cell)) %>% spread(embryo,nCell)



data.ob.umap <-  data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,embryo:sub_pred_EML,seurat_clusters)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))%>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_20),by="cell")

data.ob.umap <- data.ob.umap %>% mutate(EML=recode(seurat_clusters,'C0'="TLC","C1"="TLC","C2"="ELC"))

HB.coseq.smt2.meta.filter <- meta.filter %>% left_join(data.ob.umap %>% select(cell,EML),by="cell")
saveRDS(HB.coseq.smt2.meta.filter,paste0("tmp_data/",TD,"/HB_coseq_smt2.updated.meta.filter.rds"))
saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/HB_coseq_smt2.updated.data.ob.umap.rds"))
saveRDS(data.ob,paste0("tmp_data/",TD,"/HB_coseq_smt2.data.ob.rds"))
data.ob.umap %>% group_by(batch,embryo,EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(embryo,nCell)
