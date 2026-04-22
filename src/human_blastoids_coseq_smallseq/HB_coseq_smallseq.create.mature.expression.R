#' ---
#' title: "loading ft from each cells (HB coseq small-seq)"
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

mt.perc.cutoff <- 0.25
Total_reads.cutoff <- 0.5*10^6
miRNA.nExpG.cutoff <- 75
tRNA.nExpG.cutoff <- 50

sel.type <- c("miRNA","snoRNA","snRNA","tRNA","rRNA","piRNA","os_piRNA")

# ##Loading transcript annotation
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)


td="HB_smallseq"



FunBaiscload <- function(cn,dn) {
  temp_file=paste0("tmp_data/",dn,"/",cn,"/aligned.detailed.od.exp.txt")
  if (!file.exists(temp_file)) {
    return(NULL)
  }else if (file.info(temp_file)$size == 0) {
    return(NULL)
  }else{
    temp <-  read.delim(temp_file,stringsAsFactors = F,head=F) %>% tbl_df() %>% separate(V1,c("ID","type","len","refbp","fabp","combp","IS5p","IS3p"),sep=":") %>% rename(umi=V2) %>% mutate(type=gsub("os-piRNA","os_piRNA",type)) %>% filter(ID %in%  chr.smallRNA.id)
    temp.mature.exp <- temp %>% filter(type %in% c("miRNA","piRNA","os_piRNA")) %>% filter(len < 40) %>% group_by(ID) %>% summarise(umi=sum(umi)) %>% bind_rows( temp %>% filter(!type %in% c("miRNA","piRNA","os_piRNA"))  %>% group_by(ID) %>% summarise(umi=sum(umi)))%>% inner_join(trans.anno$mature,by="ID") %>% mutate(cell=cn) %>% ungroup()
    temp.RL.dis <- temp %>% group_by(type,len) %>% summarise(umi=sum(umi)) %>% mutate(type_len=paste0(len,"nt_",type))%>% mutate(cell=cn)%>% ungroup()
    temp.qc <- temp.mature.exp %>% group_by(type) %>% summarise(n=n()) %>% mutate(type=paste0(type,".nExpG")) %>% bind_rows(temp.mature.exp %>% group_by(type) %>% summarise(n=sum(umi))%>% mutate(type=paste0(type,".UMI")))  %>%mutate(cell=cn) %>% spread(type,n)
    temp.qc$SmallRNA_UMI.anno <- sum(temp.mature.exp$umi); temp.qc$SmallRNA_UMI.total.anno <- sum(temp$umi)
    temp.out.list <- list()
    temp.out.list$mature.exp <- temp.mature.exp%>% ungroup()
    #temp.out.list$RL.dis <- temp.RL.dis%>% ungroup()
    temp.out.list$qc <- temp.qc%>% ungroup()
    return(temp.out.list)
  }
}


meta.detail.out <- readRDS(paste0("tmp_data/HB_smallseq/HB.small.raw.meta.rds"))
table(duplicated(meta.detail.out $cell))

raw_fastq_file <- data.frame(cell=read.delim("bin/human_blastoids_batch/Sample.name.ID",head=F)$V1,batch="HB.smallseq") %>% bind_rows(data.frame(cell=read.delim("bin/HB_ME_July_2025/extra.HB_smallseq/extra.HB_smallseq.Sample.ID",head=F)$V1,batch="HB.coseq.small")) 

raw_fastq_file %>% filter(!cell %in% meta.detail.out$cell)  %>% pull(batch) %>% table()
meta.detail.out %>% filter(!cell %in% raw_fastq_file$cell)  %>% pull(batch) %>% table()

Q <- read.csv(paste0("tmp_data/",td,"/HB.small.QC.counts.csv"),row.names = 1,stringsAsFactors = F)
table(colnames(Q) %in% meta.detail.out$cell)
meta <- meta.detail.out 

temp.list <- list()
for (n in (meta %>% pull(cell))) {
  temp.list[[n]] <- n
}

check=FALSE
if (check) {
  total.input.list <- foreach (x=temp.list,n=names(temp.list),.combine=c) %dopar% {
    rv=list()
    rv[[n]]=FunBaiscload(x,"HB_smallseq" ) ## paralelly running #mast
    rv
  }
  
  counts <- list()
  counts$mature <- total.input.list %>% lapply(function(x) {x$mature.exp %>% select(-type)}) %>% do.call("bind_rows",.)%>% spread(cell,umi) %>% replace(.,is.na(.),0) 
    meta.all <- total.input.list %>% lapply(function(x) {x$qc}) %>% do.call("bind_rows",.)%>% replace(.,is.na(.),0) %>% inner_join(meta,by="cell")
  
  #' combine mapping information
  #Q <- Q %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% rename(MappedReads=MP_reads)%>% rename(AfterCutAdp_reads=beforemapping,SmallRNA_UMI=small_UMI_total) %>% select(cell,Total_reads,AfterCutAdp_reads,MappedReads,SmallRNA_UMI)
  Q <- Q %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% rename(AfterCutAdp_reads=beforemapping,MappedReads=MP_reads,SmallRNA_UMI=small_UMI_total) %>% select(cell,Total_reads,AfterCutAdp_reads,MappedReads,SmallRNA_UMI)
  chr.dis.out <- read.delim("tmp_data/HB_smallseq/HB.small.dedup.uniqMap.chr.stat.total",stringsAsFactors = F,head=F) %>% tbl_df() %>% rename(chr=V1,ct=V2,cell=V3) %>% filter(chr %in% paste0("chr",c(1:22,"X","Y","MT"))) %>% spread(cell,ct) %>% replace(.,is.na(.),0) %>% tibble::column_to_rownames("chr")
  Q <- Q %>% inner_join(as.matrix(chr.dis.out["chrMT",]/colSums(chr.dis.out)) %>% t() %>% as.data.frame()  %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% rename(mt.perc=chrMT)
  
  meta.all <- meta.all %>% left_join(Q,by="cell")
  save(counts,meta.all,file=paste0("tmp_data/",td,"/HB_coseq_smallseq.all.raw.data.Rdata"))
  saveRDS(meta.all,file=paste0("tmp_data/",td,"/HB_coseq_smallseq.all.only.meta.all.rds"))
}






