#' ---
#' title: DE miRNAs and snoRNA
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
  library(slingshot)
})


set.seed(1)
DIR <- paste0(base_dir,"/My_project/mouse_smallseq_preimp")
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)


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


#' loading gene meta data
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


sel.exp.list <- list()
sel.exp.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".norm.list.rds"))
sel.exp.list$miRNA <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

#setting
#sel.type <- c("miRNA","snoRNA")
sel.type <- c("miRNA","snoRNA","tRNA","snRNA","piRNA","rRNA")
zs.limit <- 2

#' loading data
cells.list <- meta.filter %>% filter(batch %in% c("batch1","batch2","batch4")) %>% mutate(SID=sub_EML) %>% split(.$SID) %>% lapply(function(x){x$cell})
cells.list$ICMs <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="ICM") %>% pull(cell)
cells.list$TEs <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="TE") %>% pull(cell)
cells.list$L2C <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="L2and4C") %>% filter(devTime=="2C") %>% pull(cell)
cells.list$L4C <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="L2and4C") %>% filter(devTime=="4C") %>% pull(cell)
cells.list$L8C <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="L8CM") %>% filter(devTime=="8C") %>% pull(cell)
cells.list$L16C <- meta.filter%>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML=="L8CM") %>% filter(devTime=="16C") %>% pull(cell)


lapply(cells.list,length)

#' DEG results
raw.DEG.stat=matrix(nrow=2,ncol=24)
rownames(raw.DEG.stat) <- c("up_regulated","down_regulated")
colnames(raw.DEG.stat) <- c(

  paste("sperm","vs","oocyte",sep="_"),
  paste("L2C","vs","oocyte",sep="_"),
  paste("L2C","vs","sperm",sep="_"),
  paste("L4C","vs","L2C",sep="_"),
  paste("L8C","vs","L4C",sep="_"),
  paste("L16C","vs","L8C",sep="_"),
 
  paste("Early_ICM","vs","L16C",sep="_"),
  paste("EB_TE","vs","L16C",sep="_"),
  
  paste("L8C","vs","L2C",sep="_"),
  paste("L16C","vs","L2C",sep="_"),
  paste("ICMs","vs","L2C",sep="_"),
  paste("TEs","vs","L2C",sep="_"),
  
  paste("ICMs","vs","L8C",sep="_"),
  paste("TEs","vs","L8C",sep="_"),
  
  paste("ICMs","vs","L16C",sep="_"),
  paste("TEs","vs","L16C",sep="_"),
  
  
  paste("Late_ICM","vs","Early_ICM",sep="_"),
  
  paste("MB_TE","vs","EB_TE",sep="_"),
  paste("LB_TE","vs","MB_TE",sep="_"),
  paste("LB_TE","vs","EB_TE",sep="_"),
  
  paste("Early_ICM","vs","EB_TE",sep="_"),
  paste("Early_ICM","vs","MB_TE",sep="_"),
  paste("Late_ICM","vs","LB_TE",sep="_"),
  
  paste("ICMs","vs","TEs",sep="_")
)

savefile <- paste0("tmp_data/",TD,"/lineage.segregation.DEG.out.Rdata")
DEG.results.list <- list()
DEG.stat.list <- list()


if (file.exists(savefile)){
  load(savefile,verbose = T)
}else{
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
  #' fix the "_" issues
  trans.anno$mature <- trans.anno$mature %>% mutate(ID=gsub("_","-",ID))
  trans.anno$prec <- trans.anno$prec %>% mutate(ID=gsub("_","-",ID))
  counts.filter$mature <- counts.filter$mature %>% mutate(ID=gsub("_","-",ID))
  chr.smallRNA.id <- gsub("_","-",chr.smallRNA.id)
  for (st in sel.type) {
    DEG.stat <- raw.DEG.stat
    small.ID <- trans.anno$mature %>% filter(type==st) %>% pull(ID) %>% intersect(counts.filter$mature$ID)
    chr.small.ID <- trans.anno$mature %>% filter(type==st) %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
    counts.small.filter <- (counts.filter$mature %>% filter(ID %in% chr.small.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]
    
    temp.M <- meta.filter %>% filter(batch %in% c("batch1","batch2","batch4")) %>% filter(EML!="unknown")
    temp.cells <- temp.M$cell
    temp.counts <- counts.small.filter [,temp.M$cell]
    temp.sel.expG <-  rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >=1) >=2] %>% intersect(rownames(sel.exp.list[[st]]))
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
   
    temp.norm <- sel.exp.list[[st]]
    data.ob@assays$RNA@layers$data <- as.matrix(temp.norm[temp.sel.expG,colnames(data.ob)])
    
    
    DEG.results <- list()
    for (n in colnames(DEG.stat)) {
      DEG.results[[n]] <- n
    }
    
    for (n in names(DEG.results)) {
      print(n)
      DEG.results[[n]] <- FunDEG_custom(data.ob,cells.list,n)
    }
    
    for (n in names(DEG.results)) {
      print(n)
      DEG.stat["up_regulated",n] <- nrow(DEG.results[[n]][["DEG.result.up"]])
      DEG.stat["down_regulated",n]  <- nrow(DEG.results[[n]][["DEG.result.down"]])
      print(paste(DEG.stat["up_regulated",n],DEG.stat["down_regulated",n]))
    }
    DEG.results.list[[st]] <- DEG.results
    DEG.stat.list[[st]] <- DEG.stat
  }
  save(DEG.results.list,DEG.stat.list,sel.exp.list,file=savefile)
}


