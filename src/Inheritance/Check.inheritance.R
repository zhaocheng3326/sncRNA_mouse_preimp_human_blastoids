#' ---
#' title: Inheritance
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
  library(ggtern)
  library(ggrepel)
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


# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=ifelse(EML=="prelineage" & devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(EML=="prelineage" & !devTime %in% c("2C","4C"),"L8CM",EML))  %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

C2C12MC.ov <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/C2MC_C12MC_ov_small_anno.bed",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V4 %in% c("C2MC","C12MC"))  %>% select(V4,V11,V12) %>% unique() %>% rename(cluster=V4,ID=V11,ncType=V12)
#' update the full cell annotation
meta.filter <- meta.filter  %>% left_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell") ### need to be left_join

#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap,by="cell")
meta.filter <- meta.filter %>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))


#' loading average expression
data.ob <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))
data.ob.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".data.ob.list.rds"))
data.ob.list$miRNA <- data.ob

#' raw counts
#load(paste0("tmp_data/","merge_batch_soft_link","/all.raw.data.Rdata"),verbose=T)
#rm(meta.all)
load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)
small.anno.bed <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/small.anno.bed",sep="\t",stringsAsFactors=F,header = F) %>% tbl_df()
counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))

#' loading DEG results
load(paste0("tmp_data/",TD,"/lineage.segregation.DEG.out.Rdata"),verbose=T)

snRNA.type <- c("miRNA","snoRNA","tRNA","piRNA","snRNA","rRNA")#,"Mt_tRNA"
sncRNA.list <- small.anno.bed  %>% filter(V8 %in% snRNA.type ) %>% split(.$V8) %>% lapply(function(x){x$V7 %>% unique()})
# 
for (n in snRNA.type ) {
  if (n!="Mt_tRNA") {
    sncRNA.list [[n]] <-   sncRNA.list [[n]] %>% intersect(chr.smallRNA.id)
  }
}
lapply(sncRNA.list,length)

plot.results <- list()

temp.pt.col <- c(setNames(lineage.col.set[c("sperm","oocyte","L2and4C")],c("sperm_spec","oocyte_spec","L2C_spec")),"os_spec"="royalblue","s2_spec"="#00B9E3","o2_spec"="#FB61D7","shared"="grey44","uncertained"="grey77") 


temp.col.set <- devTime.col.set[c("2C","oocyte","sperm")] %>% setNames(c("L2C","oocyte","sperm"))
temp.plot <- list()
temp.sncRNA.int.list <- list()
temp.sncRNA.deg.int.list <- list()
temp.sncRNA.nExp.detail <- list()

for (n in snRNA.type) {
  hc <- 0.25
  data.temp <-   subset(data.ob.list[[n]],cell=(data.all.ob.umap %>% filter(batch!="Split") %>% pull(cell)))#& devTime %in% c("oocyte","sperm","2C")
  temp.meta.filter <- meta.filter %>% filter(batch!="Split"  & devTime %in% c("oocyte","sperm","2C"))  #%>% filter(cell %in% rownames(data.temp@meta.data))
  temp.meta.filter$batch %>% table()
  temp.meta.filter$devTime %>% table()
  temp.exp.nC <- counts.filter$mature %>% filter(ID %in% sncRNA.list[[n]]) %>% select(ID,temp.meta.filter$cell) %>% gather(cell,ct,-ID) %>% mutate(ct=ifelse(ct >0,1,0)) %>% filter(ct >0) %>% inner_join(temp.meta.filter %>% select(cell,devTime) %>% unique() %>% group_by(devTime) %>% mutate(nBGC=n_distinct(cell)),by="cell") %>% mutate(devTime=ifelse(!devTime %in% c("oocyte","sperm"),paste0("L",devTime),devTime)) %>% group_by(ID,devTime,nBGC) %>% summarise(nC=n_distinct(cell)) %>% mutate(pct=nC/nBGC) 
  temp.plot[[n]] <- temp.exp.nC %>% split(.,.$devTime) %>% lapply(function(x){x$ID})%>% ggvenn::ggvenn(set_name_size=3,text_size = 3,show_percentage = F,fill_color = as.vector(temp.col.set ))+ggtitle(n)+FunTitle()
  temp.sncRNA.nExp.detail[[n]] <- temp.exp.nC %>% ungroup() %>% mutate(lineage=recode(devTime,"L2C"="2C"),ncType=n) %>% select(ID,lineage,ncType) %>% unique()
  
  temp.list <- list()
  temp.list$shared <- temp.exp.nC %>% group_by(ID) %>% summarise(nL=n_distinct(devTime)) %>% filter(nL==3) %>% pull(ID) %>% unique()
  temp.list$sperm_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="sperm") %>% pull(ID) %>%setdiff(temp.exp.nC %>% filter(devTime %in% c("oocyte","L2C")) %>% pull(ID) %>% unique())
  temp.list$oocyte_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="oocyte") %>% pull(ID) %>%setdiff(temp.exp.nC %>% filter(devTime %in% c("sperm","L2C")) %>% pull(ID) %>% unique())
  temp.list$L2C_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="L2C") %>% pull(ID) %>%setdiff(temp.exp.nC %>% filter(devTime %in% c("sperm","oocyte")) %>% pull(ID) %>% unique())
  temp.list$s2_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="sperm") %>% pull(ID) %>% intersect(temp.exp.nC %>% filter(pct>hc & devTime=="L2C") %>% pull(ID)) %>% setdiff(temp.exp.nC %>% filter(devTime %in% c("oocyte")) %>% pull(ID) %>% unique())
  temp.list$o2_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="oocyte") %>% pull(ID) %>% intersect(temp.exp.nC %>% filter(pct>hc & devTime=="L2C") %>% pull(ID)) %>% setdiff(temp.exp.nC %>% filter(devTime %in% c("sperm")) %>% pull(ID) %>% unique())
  temp.list$os_spec <- temp.exp.nC %>% filter(pct>hc & devTime=="oocyte") %>% pull(ID) %>% intersect(temp.exp.nC %>% filter(pct>hc & devTime=="sperm") %>% pull(ID)) %>% setdiff(temp.exp.nC %>% filter(devTime %in% c("L2C")) %>% pull(ID) %>% unique())
  temp.list$uncertained <- temp.exp.nC$ID %>% setdiff(unlist(temp.list)) %>% unique()
  for (s in names(temp.list)) {
    if (length( temp.list[[s]])!=0) {
      temp.list[[s]] <- data.frame(ID= temp.list[[s]],type=s,ncType=n) %>% tbl_df()
    }else{
      temp.list[[s]] <- NULL
    }
  }
  
  temp.sncRNA.int.list[[n]] <- temp.list %>% do.call("bind_rows",.)
  
  table(temp.sncRNA.int.list[[n]]$type)
  
  DEG.results <- DEG.results.list[[n]]
  temp.DEG <- c(DEG.results$sperm_vs_oocyte$DEG.result$gene,DEG.results$L2C_vs_oocyte$DEG.result$gene,DEG.results$L2C_vs_sperm$DEG.result$gene) %>% unique()

}
cowplot::plot_grid(plotlist=temp.plot)
plot.results$sncRNA.inherit.raw.venn <- temp.plot

scnRNA.inh.out <- temp.sncRNA.int.list %>% do.call("bind_rows",.)
scnRNA.inh.out %>% group_by(ncType,type) %>% summarise(nID=n_distinct(ID)) %>% spread(type,nID)

#' stat
temp.stat <- scnRNA.inh.out %>% group_by(type,ncType) %>% summarise(nGene=n()) %>% group_by(ncType) %>% mutate(prop=nGene/sum(nGene)) %>% select(-nGene)  %>% spread(ncType,prop)%>% replace(.,is.na(.),0) %>% gather(ncType,prop,-type) %>% mutate(prop=100*prop)
plot.results$sncRNA.inherit.prop <- temp.stat %>% mutate(ncType=factor(ncType,rev(snRNA.type),ordered = T)) %>% mutate(type=factor(type,rev(names(temp.pt.col)),ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=ncType,y=prop,fill=type),stat="identity")+coord_polar(theta = "y") +scale_fill_manual(values=temp.pt.col)+theme_classic()
plot.results$sncRNA.inherit.prop 


