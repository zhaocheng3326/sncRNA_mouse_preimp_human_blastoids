#' ---
#' title: WGCNA for coseq mRNA
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
  #condaENV <- "/home/chenzh/miniconda3/envs/R4.3" 
  #print(condaENV)
  #LBpath <- paste0(condaENV ,"/lib/R/library")
  #.libPaths(LBpath)
  base_dir="/home/chenzh"
}



suppressMessages({
  library(Seurat)
  library(scran)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  library(WGCNA)
  library(hdWGCNA)
  library(EnsDb.Mmusculus.v79)
  
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
MSplit.meta.filter <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds"))

wgcna_name <- "coseq_smt2"


#' co-smt2 dataset only
if (file.exists(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))) {
  data.ref <- readRDS(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))
  data.ref.perum <- readRDS(paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ref.perum.ob.rds"))
}else{
  temp.M <- MSplit.meta.filter %>% filter(pj=="co_smt2")
  temp.sel.expG <- rownames(MSplit.counts.filter)[rowSums(MSplit.counts.filter[,c(temp.M$cell)] >=1) >2]
  temp.counts <- MSplit.counts.filter[temp.sel.expG,temp.M$cell]
  
  nGene=2000;npc=25; data.ob <- CreateSeuratObject(temp.counts, meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  #temp.sce <-  SingleCellExperiment(list(counts=as.matrix(temp.counts)),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
  #temp.norm <- scuttle::normalizeCounts(temp.sce)
  #data.ob@assays$RNA@layers$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.ob)])
  data.ob <- data.ob %>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:npc,verbose=F) #%>% FindNeighbors( dims = 1:npc,verbose = FALSE) %>%  FindClusters(resolution = 0.8,verbose = FALSE)  #,vars.to.regress=c("nGene") 
  
  data.temp <- subset(data.ob,cells=(temp.M %>% filter(EML!="unknown") %>% pull(cell)))
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("EML"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("pj"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("stage"),label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,group.by=c("devTime"),label=T)+NoAxes()+NoLegend()
    #FeaturePlot(data.temp,c("CCR7"))+NoAxes()+NoLegend()
  )
  
  
  #' setup for hdWGCNA
  data.temp <- SetupForWGCNA(data.temp,gene_select = "fraction", fraction =  round(3/nrow(data.temp@meta.data),2), wgcna_name = wgcna_name )
  
  #' group cells
  data.temp <- MetacellsByGroups(seurat_obj = data.temp, group.by = c("EML"),  reduction = 'pca', k = 10, max_shared = 5, ident.group = 'EML',min_cells=10)
  data.temp <- NormalizeMetacells(data.temp)
  
  # using scran for merged object
  #data.merge.temp <- GetMetacellObject(data.temp, wgcna_name)
  #temp.sce <-  SingleCellExperiment(list(counts=data.merge.temp@assays$RNA$counts),colData=data.merge.temp@meta.data) %>% computeSumFactors()
  #temp.norm <- scuttle::normalizeCounts(temp.sce)
  #data.merge.temp@assays$RNA@layers$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.merge.temp)])
  #data.temp <- SetMetacellObject(data.temp, data.merge.temp, wgcna_name)
  
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
  data.temp <- ConstructNetwork(data.temp,tom_outdir=paste0('tmp_data/',TD,'/TOM/'),overwrite_tom = TRUE,tom_name =wgcna_name)
  
  #' compute all MEs, cell x Eigengene
  data.temp <- ModuleEigengenes(data.temp)
  
  #' compute eigengene-based connectivity (kME):
  data.temp <- ModuleConnectivity(data.temp)
  
  data.ref <- data.temp
  saveRDS(data.ref,paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ob.rds"))
  
  
  
  md.list <- list()
  md.list$MEs <- GetMEs(data.ref) 
  md.list$modules <- GetModules(data.ref) %>% subset(module != 'grey')
  md.list$md.anno <- md.list$modules %>% select(module) %>% unique() %>% tbl_df() %>% tibble::rowid_to_column("gene_md")  %>% mutate(gene_md=paste("gene_md",gene_md,sep="_")) %>% rename(gene_mdc=module)
  
  saveRDS(md.list,paste0("tmp_data/",TD,"/TOM/",wgcna_name,".md.list.rds"))
  
  #' permutation calcuate the Z-score
  set.seed(123)
  data.ref.perum <- data.ref
  data.ref.perum <- ProjectModules(data.ref.perum,data.temp,wgcna_name = wgcna_name,wgcna_name_proj="perum", assay="RNA" )
  data.ref.perum <- ModulePreservation(SetDatExpr(data.ref.perum),SetDatExpr(data.ref),name="perum",verbose=3,n_permutations=200)
  saveRDS(data.ref.perum,paste0("tmp_data/",TD,"/TOM/",wgcna_name,".data.ref.perum.ob.rds"))
  
}

