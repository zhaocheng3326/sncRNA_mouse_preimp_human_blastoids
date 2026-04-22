#' ---
#' title: check the small RNA part for split seq (based on miRNA expression)
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
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds")) 
counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]

co.mk <- c("mmu-miR-99a-5p","mmu-miR-320-3p","mmu-miR-27a-3p","mmu-let-7b-5p","mmu-miR-22-3p","mmu-miR-23b-3p","mmu-miR-23a-3p","mmu-miR-221-3p", "mmu-miR-205-3p","mmu-miR-203-3p","mmu-miR-144-3p","mmu-miR-181a-5p","mmu-miR-205-5p","mmu-miR-130a-3p")


#' split batch rename
meta.filter <- meta.filter %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split")) %>% mutate(pj=batch)


#' running
sel.type <- c("miRNA")
zs.limit <- 2
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)


#'loading data (for coseq small only)
meta.filter <- meta.filter %>% filter(batch=="Split")
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]

#' loading coseq smt2 cell annotation
MSplit.meta.filter <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds")) 
#meta.filter <- meta.filter %>% left_join(MSplit.meta.filter %>% select(cell,EML),by="cell")  %>% mutate(EML=ifelse(is.na(EML),"None",EML))


temp.M <- meta.filter %>% filter(batch=="Split")%>%   left_join(readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.meta.filter.mod.withsex.rds")) %>% select(embryo,sex) %>% unique(),by="embryo")
temp.counts <- counts.miRNA.filter[,temp.M$cell]
temp.M <- temp.M %>% left_join(MSplit.meta.filter %>% select(cell,EML),by="cell")  %>% mutate(EML=ifelse(is.na(EML),"None",EML))
temp.cells <- temp.M$cell
temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2]

if (file.exists(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds"))
  mk.list <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.split.mk.list.rds"))
  fm.list <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.split.fm.list.rds"))
  sel.exp <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.norm.rds"))
}else{
  
  # npc <- 25;ngene <- length(temp.sel.expG)
  npc <- 10;ngene <- 150
  data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  temp.norm <- scuttle::normalizeCounts(temp.sce)
  data.ob@assays$RNA@layers$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.ob)])
  data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE,return.model = TRUE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
  # data.ob <- data.ob %>% FindVariableFeatures(verbose=F,nfeatures=ngene) %>% ScaleData(verbose = FALSE) %>% RunPCA( verbose = FALSE) %>% RunUMAP( reduction = "pca", dims = 1:npc,verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:npc,verbose = FALSE,k.param=10)
  
  
  data.temp <- data.ob  %>% FindClusters( resolution = 0.3,verbose = FALSE)
  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "sex",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "stage",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% mutate(RNA_EML=EML)%>% select(c(cell,embryo,sex,devTime,stage,RNA_EML)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp)))) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_15),by="cell")
  
  #data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% mutate(RNA_EML=EML)%>% select(c(cell,Total_reads:SmallRNA_UMI.anno,RNA_EML)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp)))) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_15),by="cell")
  data.ob.umap <- data.ob.umap %>% mutate(small_EML=recode(seurat_clusters,"C1"="prelineage","C3"="prelineage","C0"="TE","C2"="ICM"))  %>% mutate(EML=ifelse(RNA_EML==small_EML,RNA_EML,"unknown"))%>% mutate(EML=ifelse(RNA_EML %in% c("None","unknown"),small_EML,EML)) %>% mutate(EML=ifelse(devTime %in% c("64C") & EML=="prelineage","unknown",EML))%>% mutate(EML=ifelse(devTime %in% c("8C","4C","2C") & EML %in% c("ICM","TE"),"unknown",EML)) #%>% mutate(small_EML=ifelse(devTime %in% c("8C","16C"),"prelineage",small_EML))
  
  
  #Idents(data.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"seurat_clusters"])
  data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
  data.ob@meta.data$small_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"small_EML"]
  data.ob@meta.data$RNA_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"RNA_EML"]
  data.ob@meta.data$seurat_clusters <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"seurat_clusters"]
  
  #' get the markers for different clusters
  mk.list <- list()
  fm.list <- list()
  temp.cells <- data.ob.umap %>% filter(EML %in% c("prelineage","ICM","TE")) %>% pull(cell)
  data.deg <- subset(data.ob,cells=temp.cells)
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.deg),"EML"])
  mk.list$EML_mk <- FunRF_FindAllMarkers_para(data.deg)
  fm.list$EML_mk <- FindAllMarkers(data.deg,only.pos = T) %>% tbl_df()
  
  data.deg <- subset(data.ob,cells=temp.cells)
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.deg),"EML"])
  data.deg <- RenameIdents(data.deg,"ICM"="ICM_TE","TE"="ICM_TE")
  mk.list$EML_mk$sig <-   mk.list$EML_mk$sig %>% bind_rows((FunRF_FindAllMarkers_para(data.deg))$sig %>% tbl_df() %>% filter(set=="ICM_TE"))
  
  fm.list$EML_mk <-   fm.list$EML_mk %>% bind_rows(FindAllMarkers(data.deg,only.pos = T) %>% tbl_df() %>% filter(cluster=="ICM_TE"))
  
  temp.cells <- data.ob.umap  %>% pull(cell) 
  data.deg <- subset(data.ob,cells=temp.cells)
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"seurat_clusters"])
  fm.list$cluster_mk <-   FindAllMarkers(data.deg,only.pos = T) %>% tbl_df() 
  
  
  temp.cells <- data.ob.umap  %>% pull(cell) 
  data.deg <- subset(data.ob,cells=temp.cells)
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"seurat_clusters"])
  #data.deg <- RenameIdents(data.deg,"C2"="ICM","C0"="TE","C1"="TE")
  mk.list$cluster_mk <- FunRF_FindAllMarkers_para(data.deg)
  
  mk.list$EML_mk$sig %>% full_join(mk.list$cluster_mk$sig,by=c("gene","set")) %>% arrange(desc(power.y))
  mk.list$EML_mk$sig %>% group_by(set) %>% top_n(20,power)%>% inner_join(mk.list$cluster_mk$sig %>% group_by(set) %>% top_n(20,power) ,by=c("gene")) %>% arrange(desc(power.y)) %>% pull(set.x) %>% table()
  
  
  meta.filter.updataed <- meta.filter %>% left_join(data.ob.umap %>% select(cell,EML,RNA_EML,small_EML),by="cell")
  
  
  temp.cells <- data.ob.umap %>% filter(EML %in% c("ICM","TE")) %>% pull(cell)
  data.deg <- subset(data.ob,cells=temp.cells)
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.deg),"EML"])
  temp.DEG <- FindMarkers(data.deg,ident.1="ICM",ident.2="TE") %>% tibble::rownames_to_column("gene")%>% tbl_df() %>% mutate(fdr=p.adjust(p_val,method="BH",n=nrow(data.deg)))
  saveRDS(temp.DEG,paste0("tmp_data/",TD,"/mouse_coseq_smallseq.updated.ICMvsTE.DEG.rds"))

  
  #' save object
  saveRDS(data.ob,paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds"))
  saveRDS(meta.filter.updataed,paste0("tmp_data/",TD,"/Msmall.coseq_small.updated.meta.rds"))
  saveRDS(mk.list,paste0("tmp_data/",TD,"/Msmall.coseq_small.split.mk.list.rds"))
  saveRDS(fm.list,paste0("tmp_data/",TD,"/Msmall.coseq_small.split.fm.list.rds"))
  saveRDS(temp.norm,paste0("tmp_data/",TD,"/Msmall.coseq_small.norm.rds"))
}
