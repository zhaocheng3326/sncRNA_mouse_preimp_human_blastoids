#' ---
#' title: strict DE detection among blastocyst/blastoids/naive primed cells
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


raw.DEG.stat=matrix(nrow=2,ncol=7)
rownames(raw.DEG.stat) <- c("up_regulated","down_regulated")
colnames(raw.DEG.stat) <- c(
  
  paste("ICM","vs","TE",sep="_"),
  paste("ELC","vs","TLC",sep="_"),
  paste("naive","vs","primed",sep="_"),
  
  paste("ICM","vs","ELC",sep="_"),
  paste("TE","vs","TLC",sep="_"),
  
  paste("ICM","vs","naive",sep="_"),
  paste("ELC","vs","naive",sep="_")
)

savefile <- paste0("tmp_data/",TD,"/HB.HE.NP.DEG.out.strict.Rdata")
DEG.results.list <- list()
DEG.stat.list <- list()
sel.exp.list <- list()

sel.type <- c("miRNA","snoRNA","tRNA","snRNA","piRNA","rRNA")



if (file.exists(savefile)){
  load(savefile,verbose = T)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB.HE.NP.DEG.meta.filter.strict.rds"))
}else{
  #' human small meta
  load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)
  
  
  HB.meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
  HB.counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))
  
  
  HE.meta.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.coord.check.rds")%>% mutate(rename_EML=ifelse(EML %in% c("EarlyTE","mural","polar"),"TE",EML)) %>% mutate(EML=rename_EML) %>% select(cell,embryo,devTime,batch,EML)
  HE.counts.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")
  
  NP.meta.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/raw_data/GSE81287/small.meta.filter.rds") %>% mutate(embryo=devTime) %>% select(cell,embryo,devTime,batch,EML) %>% filter(EML!="HEK")
  NP.counts.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/raw_data/GSE81287/small.counts.filter.rds")
  
  
  meta.filter <- HB.meta.filter %>% select(cell,embryo,batch,EML,devTime) %>% bind_rows(HE.meta.filter %>% select(cell,embryo,batch,EML,devTime) )%>% bind_rows(NP.meta.filter %>% select(cell,embryo,batch,EML,devTime) )
  table(duplicated(meta.filter$cell))
  table(meta.filter$batch)
  
  meta.filter <- meta.filter %>% filter(!batch %in% c("HB.coseq.small","HumanSplit")) 
  cells.list <- meta.filter %>% filter(devTime!="E5") %>% mutate(SID=EML) %>% split(.$SID) %>% lapply(function(x){x$cell})
  
  for ( st in sel.type) {
    DEG.stat <- raw.DEG.stat
    miRNA.ID <- trans.anno$mature %>% filter(type==st) %>% pull(ID) %>% intersect(HB.counts.filter$mature$ID)%>% intersect(HE.counts.filter$mature$ID)%>% intersect(NP.counts.filter$mature$ID)
    chr.miRNA.ID <- miRNA.ID  %>% intersect(chr.smallRNA.id) 
    counts.miRNA.filter <- HB.counts.filter$mature %>% filter(ID %in% chr.miRNA.ID)  %>% inner_join(HE.counts.filter$mature %>% filter(ID %in% chr.miRNA.ID),by="ID") %>% inner_join(NP.counts.filter$mature %>% filter(ID %in% chr.miRNA.ID),by="ID") %>% tibble::column_to_rownames("ID")
    counts.miRNA.filter <- counts.miRNA.filter[,meta.filter$cell]       
    
    rownames(counts.miRNA.filter) <- gsub("_","-",rownames(counts.miRNA.filter))
    
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
    mBN.sce.ob <- multiBatchNorm(sce.ob$SB,sce.ob$B1,sce.ob$TB,sce.ob$A1,sce.ob$HB.smallseq,sce.ob$batchNPH)
    names(mBN.sce.ob) <- c("SB","B1","TB","A1","HB.smallseq","batchNPH")
    mBN.sce.ob %>% lapply(function(x) {data.frame(cell=colnames(x),sf=sizeFactors(x)) %>% tbl_df() %>% return()})  %>% do.call("bind_rows",.) %>% inner_join(meta.filter %>% select(cell,batch)) %>% mutate(od=batch,ordered = T) %>% ggplot()+geom_violin(mapping=aes(x=od,y=sf),scale = "width",fill="royalblue")+theme_classic()
    sel.exp <- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
    sel.exp.list[[st]] <- sel.exp
    
    
    temp.M <- meta.filter %>% filter(cell %in% colnames(sel.exp ))
    temp.counts <- counts.miRNA.filter[,temp.M$cell]
    temp.sel.expG <- rownames(sel.exp)  #%>% setdiff(chrX.miRNA.unique)
    data.ob <- CreateSeuratObject(temp.counts[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    data.ob@assays$RNA$data <- as.matrix(sel.exp[temp.sel.expG,colnames(data.ob)])
    
    DEG.results <- list()
    for (n in colnames(DEG.stat)) {
      DEG.results[[n]] <- n
    }
    
    for (n in names(DEG.results)) {
      print(n)
      if (n %in% c("ICM_vs_TE","ELC_vs_TLC","naive_vs_primed")) {
        DEG.results[[n]] <- FunDEG_custom(data.ob,cells.list,n)
      }else{
        DEG.results[[n]] <- FunDEG_custom_strict(data.ob,cells.list,n)
      }
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
  saveRDS(meta.filter,paste0("tmp_data/",TD,"/HB.HE.NP.DEG.meta.filter.strict.rds"))
}



