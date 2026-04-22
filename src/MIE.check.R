#' ---
#' title: "check MIE"
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




para_cal=TRUE
if (para_cal) {
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  numCores <- 11
  registerDoParallel(numCores)
}

#' miRNA family 
miRNA.family <- readRDS(paste0("tmp_data/",TD,"/miRNA.family.rds"))
miRNA.bed <- readRDS(paste0("tmp_data/",TD,"/miRNA.bed.rds"))




# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=ifelse(EML=="prelineage" & devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(EML=="prelineage" & !devTime %in% c("2C","4C"),"L8CM",EML))  %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

#' update the full cell annotation
meta.filter <- meta.filter  %>% left_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell") ### need to be left_join

#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap,by="cell")


#' markers
miRNA.FM <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.FM.mk.rds"))
miRNA.mk.out <- readRDS(file=paste0("tmp_data/",TD,"/","allCells.miRNA.mk.rds"))
miRNA.ac.sel.exp <-  readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))
data.ob <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))

data.ob.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.list.rds"))
data.ob.list$miRNA <- data.ob


plot.results <- list()
library(infotheo)
temp.plot <- list()
for (n in c("miRNA","snoRNA","snRNA","tRNA","rRNA","piRNA")) {
  data.temp <- subset(data.ob.list[[n]],cell=(data.all.ob.umap %>% filter(batch!="Split") %>% pull(cell)))
  temp_ave_exp <- AverageExpression(data.temp,group.by = "devTime")$RNA %>% as.data.frame() 
  #temp_ave_exp <- temp_ave_exp[temp_ave_exp %>% apply(1,max) >1,]
  temp_n <- ncol(temp_ave_exp[,])
  temp_mi_matrix <- matrix(0, nrow = temp_n, ncol = temp_n)
  
  for (i in 1:temp_n) {
    for (j in i:temp_n) {
      temp_mi <- mutinformation(round(temp_ave_exp [,i]), round(temp_ave_exp [,j]))
      temp_mi_matrix[i, j] <- temp_mi
      temp_mi_matrix[j, i] <- temp_mi
    }
  }
  colnames(temp_mi_matrix) <- colnames(temp_ave_exp )
  rownames(temp_mi_matrix) <- colnames(temp_ave_exp )
  temp.plot[[n]] <- temp_mi_matrix %>% pheatmap::pheatmap(scale="none",main=n,col=rev(RColorBrewer::brewer.pal(n =11, name = "RdBu"))[c(6:11)],width = 4.5,height=4.5)%>%ggplotify::as.ggplot()
}
cowplot::plot_grid(plotlist = temp.plot)
plot.results$stage.MIE.ph <- temp.plot
"tmp_data/Fig_pdf/temp.stage.MIE.check.pdf"


