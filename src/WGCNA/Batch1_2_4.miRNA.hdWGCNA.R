#' ---
#' title: WGCNA for batch1,2,4 miRNA and project for co-seq
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
  library(igraph)
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
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
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


theme_set(theme_cowplot())
wgcna_name <- "all_small_nonSplit_miRNA"



# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=recode(EML,"prelineage"="L8CM")) %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the full cell annotation
meta.filter <- meta.filter  %>% left_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap,by="cell")



data.ob <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))
sel.norm <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

#mmu-miR-34c-3p
#'
if (file.exists(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))) {
  data.ref <- readRDS(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))
  data.ref.perum <- readRDS(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ref.perum.ob.rds"))
  data.query.perum <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.ob.rds"))
}else{
  #' exclude split-seq cells
  data.sub.ob <- subset(data.ob,cell=meta.filter %>% filter(! batch %in% c("Split1","Split2") & EML!="unknown") %>% pull(cell))
  data.sub.ob@meta.data$sub_EML <- (meta.filter %>%tibble::column_to_rownames("cell"))[rownames(data.sub.ob@meta.data),"sub_EML"]
  table(data.sub.ob@meta.data$sub_EML)
  
  data.temp <- data.sub.ob
  cowplot::plot_grid(
    DimPlot(data.temp,group.by=c("sub_EML"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("EML"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("pj"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("stage"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("devTime"),label=T)+NoAxes()+NoLegend()
    #FeaturePlot(data.temp,c("CCR7"))+NoAxes()+NoLegend()
  )
  
  #' setup for hdWGCNA
  data.temp <- SetupForWGCNA(data.temp,gene_select = "fraction", fraction = round(10/nrow(data.temp@meta.data),2), wgcna_name = wgcna_name )
  #' group cells
  data.temp <- MetacellsByGroups(seurat_obj = data.temp, group.by = c("EML"),  reduction = 'pca', k = 10, max_shared = 5, ident.group = 'EML',min_cells=10)#target_metacells=ceiling(nrow(data.temp@meta.data)/10)
  #data.temp <- NormalizeMetacells(data.temp)
  # 'using scran for merged object
  data.merge.temp <- GetMetacellObject(data.temp, wgcna_name) %>% NormalizeData(verbose = F)
  table( data.merge.temp@meta.data$EML)
  temp.sce <-  SingleCellExperiment(list(counts=data.merge.temp@assays$RNA$counts),colData=data.merge.temp@meta.data) %>% computeSumFactors()
  temp.norm <- scuttle::normalizeCounts(temp.sce)
  data.merge.temp@assays$RNA@layers$data <- as.matrix(temp.norm[rownames(data.merge.temp@assays$RNA$counts),colnames(data.merge.temp)])
  data.temp <- SetMetacellObject(data.temp, data.merge.temp)
  
  #' can be used for subset of cells
  data.temp <- SetDatExpr(data.temp,assay = 'RNA', layer = 'data')
  
  #' IMP, Select soft-power threshold
  set.seed(123)
  data.temp <- TestSoftPowers(data.temp,networkType = 'signed' )
  cowplot::plot_grid(plotlist = PlotSoftPowers(data.temp ))
  
  #The general guidance for WGCNA and hdWGCNA is to pick the lowest soft power threshold that has a Scale Free Topology Model Fit greater than or equal to 0.8
  #check power
  power_table <- GetPowerTable(data.temp)
  head(power_table,20)
  
  #'construct co-expression network
  data.temp <- ConstructNetwork(data.temp,tom_outdir=paste0('tmp_data/',TD,'/TOM/'),tom_name =wgcna_name,overwrite_tom = TRUE,minModuleSize=5)
  
  #' compute all MEs, cell x Eigengene
  data.temp <- ModuleEigengenes(data.temp)
  
  #' compute eigengene-based connectivity (kME):
  data.temp <- ModuleConnectivity(data.temp)
  
  #" calculate module trait
  #data.temp@meta.data$sub_EML <- factor(data.temp@meta.data$sub_EML,c("sperm","oocyte","L2and4C","L8CM","ICM","EPI","PE","EB_TE","MB_TE","LB_TE"),ordered = T)
  #data.temp@meta.data$EML <- factor(data.temp@meta.data$EML,c("sperm","oocyte","L2and4C","L8CM","TE","ICM"),ordered = T)
  #data.temp@meta.data$devTime <- factor(data.temp@meta.data$devTime,c("sperm","oocyte","2C","4C","8C","16C","32C","64C"),ordered = T)
  data.temp@meta.data$devTime <- factor(data.temp@meta.data$devTime)
  data.temp@meta.data$EML <- factor(data.temp@meta.data$EML)
  data.temp@meta.data$sub_EML <- factor(data.temp@meta.data$sub_EML)
  
  
  data.temp  <- ModuleTraitCorrelation(data.temp,c("devTime","EML","sub_EML"))
  
  data.ref <- data.temp
  saveRDS(data.ref,paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))
  md.list <- list()
  md.list$TOM <- GetTOM(data.ref)
  md.list$MEs <- GetMEs(data.ref) 
  md.list$modules <- GetModules(data.ref) %>% subset(module != 'grey')
  md.list$mt_cor <- GetModuleTraitCorrelation(data.ref)
  md.list$md.anno <- md.list$modules %>% select(module) %>% unique() %>% tbl_df() %>% tibble::rowid_to_column("miRNA_md")  %>% mutate(miRNA_md=paste("miRNA_md",miRNA_md,sep="_")) %>% rename(miRNA_mdc=module)
  saveRDS(md.list,file=paste0("tmp_data/",TD,"/TOM/",wgcna_name,".md.list.rds"))
  
  
  #' permutation calcuate the Z-score
  set.seed(123)
  data.ref.perum <- data.ref
  data.ref.perum <- ProjectModules(data.ref.perum,data.temp,wgcna_name = wgcna_name,wgcna_name_proj="perum", assay="RNA" )
  data.ref.perum <- ModulePreservation(SetDatExpr(data.ref.perum),SetDatExpr(data.ref),name="perum",verbose=3,n_permutations=200)
  saveRDS(data.ref.perum,paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ref.perum.ob.rds"))
  
  
  #"check the split (projection)
  temp.M <- meta.filter %>% filter(batch %in% c("Split1","Split2") & EML!="unknown") %>% mutate(EML=ifelse(EML %in% c("TE","ICM"),EML,devTime)) %>% mutate(EML=recode(EML,"2C"="L2and4C","4C"="L2and4C","8C"="L8CM","16C"="L8CM","32C"="L8CM"))
  data.query <-  subset(data.ob,cell=temp.M$cell)
  data.query@meta.data$EML <- (temp.M %>%tibble::column_to_rownames("cell"))[rownames(data.query@meta.data),"EML"]
  data.query@meta.data$sub_EML <- (temp.M %>%tibble::column_to_rownames("cell"))[rownames(data.query@meta.data),"EML"]
  #data.query <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.rds"))
  #data.query <-  subset(data.query,cell=meta.filter %>% filter(batch=="Split" & small_EML!="unknown") %>% pull(cell))
  #data.query@meta.data$sub_EML <- (meta.filter %>%tibble::column_to_rownames("cell"))[rownames(data.query@meta.data),"small_EML"]
  #data.query@meta.data$sub_EML[data.query@meta.data$sub_EML=="prelineage"]="L8CM"
  #
  
  table(data.query@meta.data$sub_EML)
  data.query <- SetupForWGCNA(data.query,gene_select = "fraction", fraction = round(3/nrow(data.query@meta.data),2), wgcna_name = "coseq" )
  data.query <- MetacellsByGroups(seurat_obj = data.query, group.by = c("EML"),  reduction = 'pca', k = 10, max_shared = 5, ident.group = 'EML',min_cells=10)
  data.merge.temp <- GetMetacellObject(data.query, "coseq") %>% NormalizeData(verbose = F)
  table( data.merge.temp@meta.data$EML)
  temp.sce <-  SingleCellExperiment(list(counts=data.merge.temp@assays$RNA$counts),colData=data.merge.temp@meta.data) %>% computeSumFactors()
  temp.norm <- scuttle::normalizeCounts(temp.sce)
  data.merge.temp@assays$RNA@layers$data <- as.matrix(temp.norm[rownames(data.merge.temp@assays$RNA$counts),colnames(data.merge.temp)])
  data.query <- SetMetacellObject(data.query, data.merge.temp)
  
  set.seed(123)
  data.query <- ProjectModules( data.query,data.ref,wgcna_name = wgcna_name,wgcna_name_proj="coseq",assay="RNA" )
  setdiff(colnames(GetMEs(data.ref)),colnames(GetMEs(data.query)))
  data.query.perum <- ModulePreservation(SetDatExpr(data.query),SetDatExpr(data.ref),name="perum",verbose=3,n_permutations=200)
  saveRDS(data.query.perum,paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.ob.rds"))
  
  md.list <- list()
  md.list$MEs <- GetMEs(data.query.perum) 
  md.list$modules <- GetModules(data.query.perum) %>% subset(module != 'grey')
  md.list$md.anno <- readRDS(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".md.list.rds"))$md.anno
  saveRDS(md.list,file=paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
  
}
