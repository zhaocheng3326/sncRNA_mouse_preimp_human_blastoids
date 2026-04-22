#' ---
#' title: number of expressed miRNAs
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


#' human small meta
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T)
human.miRNA.ID <- trans.anno$mature  %>% filter(type=="miRNA")%>% pull(ID)
human.norm.exp <- readRDS(paste0("~/My_project/Extra_smncRNA/tmp_data/","Apr_2023","/","miRNA",".norm.rds"))$mBN

load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)
mouse.miRNA.ID <- trans.anno$mature  %>% filter(type=="miRNA") %>% pull(ID)
mouse.norm.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))



sel.type <- c("miRNA")

HB.meta.filter <- readRDS(paste0("tmp_data/",TD,"/HB.small.ac_miRNA.data.ob.umap.rds"))
HB.counts.filter <- readRDS(paste0("tmp_data/",TD,"/HB_coseq_smallseq.small.counts.filter.rds"))


HE.meta.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.coord.check.rds")%>% mutate(rename_EML=ifelse(EML %in% c("EarlyTE","mural","polar"),"TE",EML)) %>% mutate(EML=rename_EML) %>% select(cell,embryo,devTime,batch,EML)
HE.counts.filter <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")


ME.counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
ME.meta.filter <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update((readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=ifelse(EML=="prelineage" & devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(EML=="prelineage" & !devTime %in% c("2C","4C"),"L8CM",EML))  %>% select(cell,RNA_EML,small_EML,EML)) %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")

#meta.filter <- HB.meta.filter %>% select(cell,embryo,batch,EML,devTime) %>% bind_rows(HE.meta.filter %>% select(cell,embryo,batch,EML,devTime) ) %>% filter(!batch %in% c("HB.coseq.small","HumanSplit")) 

miRNA.sort.nExp.miRNA.list  <- list()

#'
temp.cells <-  HE.meta.filter %>% filter(devTime %in% c("E3","E4","E5","E6")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.E3E6 <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1] %>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))


temp.cells <-  HE.meta.filter %>% filter(devTime %in% c("E3","E4","E5","E6","E7")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.E3E7 <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

temp.cells <-  HE.meta.filter %>% filter(devTime %in% c("E3","E4","E5","E6","E7") & EML %in% c("ICM","TE")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.E6E7.ICM_TE <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

temp.cells <-  HE.meta.filter %>% filter(EML %in% c("ICM","TE")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.ICM_TE <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

temp.cells <-  HE.meta.filter %>% filter(EML %in% c("ICM")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.ICM <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

temp.cells <-  HE.meta.filter %>% filter(EML %in% c("TE")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- HE.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$human.TE <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

#' 
temp.cells <-  HB.meta.filter %>% filter(EML %in% c("ELC","TLC")) %>% filter(batch!="HB.smallseq") %>% pull(cell)
temp.counts <- HB.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$HB.ELC_TLC <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID) #%>% intersect(rownames(human.norm.exp))

temp.cells <-  HB.meta.filter %>% filter(EML %in% c("ELC")) %>% filter(batch!="HB.smallseq") %>% pull(cell)
temp.counts <- HB.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$HB.ELC <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID)

temp.cells <-  HB.meta.filter %>% filter(EML %in% c("TLC")) %>% filter(batch!="HB.smallseq") %>% pull(cell)
temp.counts <- HB.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$HB.TLC <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(human.miRNA.ID)

#'
temp.cells <-  ME.meta.filter %>% filter(!devTime %in% c("oocyte","sperm","2C","4C")) %>% filter(batch!="Split") %>% pull(cell)
temp.counts <- ME.counts.filter$mature %>%tibble::column_to_rownames("ID")
miRNA.sort.nExp.miRNA.list$mouse.8C_64C <- rownames(temp.counts)[rowSums(temp.counts[,temp.cells] >0) >=1]%>% intersect(mouse.miRNA.ID) #%>% intersect(rownames(mouse.norm.exp))


miRNA.sort.nExp.miRNA.list %>% lapply(head)
miRNA.sort.nExp.miRNA.list %>% lapply(length) %>% as.data.frame() %>% t() %>% as.data.frame() %>% setNames("Number of expressed miRNA")

saveRDS(miRNA.sort.nExp.miRNA.list,file=paste0("tmp_data/",TD,"/","miRNA.sort.nExp.miRNA.list.rds"))



#' plot results
plot.results <- list()
miRNA.sort.nExp.miRNA.list %>% names()
temp.plot <- list()
#temp.plot$t1 <-  miRNA.sort.nExp.miRNA.list[c("human.ICM_TE","HB.ELC_TLC")] %>% ggvenn::ggvenn(set_name_size=3,text_size = 3)
temp.plot$t2 <-  miRNA.sort.nExp.miRNA.list[c("human.ICM","HB.ELC")] %>% ggvenn::ggvenn(show_percentage = F)
temp.plot$t3 <-  miRNA.sort.nExp.miRNA.list[c("human.TE","HB.TLC")] %>% ggvenn::ggvenn(show_percentage = F)
plot.results$HE.HB.nExp.venn <- temp.plot
cowplot::plot_grid(plotlist=temp.plot,ncol=2)



#' check CMC percentage boxplot
temp.plot <- list()
#' check the CMC one mouse
temp.M <-  ME.meta.filter%>% filter(batch!="Split" & EML!="unknown") %>% mutate(SID=devTime)%>% mutate(devTime=recode(devTime,"64C"="MB","32C"="EB"))%>% mutate(SID=ifelse(EML %in% c("ICM","TE"),sub_EML,SID)) %>% mutate(SID=recode(SID,"Early_ICM"="EB_ICM","Late_ICM"="MB_ICM", "8C"="L8C","16C"="L16C","32C"="L32C","EB_TE"="Emergent_TE","MB_TE"="Early_TE","LB_TE"="Middle_TE"))  %>% filter(SID!="Emergent_TE") %>% mutate(SID=recode(SID,"Middle_TE"="MB_TE","Early_TE"="EB_TE"))
C2C12MC.ov <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/C2MC_C12MC_ov_small_anno.bed",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V4 %in% c("C2MC","C12MC")) %>% filter(V12=="miRNA") %>% select(V4,V11,V12) %>% unique()

temp.plot$mouse <- ME.counts.filter$mature %>% inner_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11),by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% right_join(temp.M %>% select(cell,SID),by="cell")  %>% select(CMC,SID,cell,ID) %>% unique()%>% group_by(SID,cell,CMC) %>% summarise(nExpMiRNA=n_distinct(ID))  %>% right_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11) %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")) %>%  mutate(SID=factor(SID,c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"),ordered = T)) %>% ggplot()+geom_boxplot(mapping=aes(x=CMC,y=nExpMiRNA/nTotalMiRNA,fill=SID),position="dodge",outlier.size = 0.25)+xlab("")+scale_fill_manual(values=c(lineage.col.set[c("L8CM","Prelineage","Early_ICM","Late_ICM","EB_TE","LB_TE")] %>% setNames(c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"))))+theme_classic()


#' check the CMC one human
temp.M <-  HE.meta.filter %>% filter(batch!="HumanSplit" & devTime!="E7") %>% mutate(SID=EML) %>% mutate(devTime=recode(devTime,"E6"="MB","E7"="LB","E5"="EB"))%>% mutate(SID=ifelse(SID %in% c("ICM","TE"),paste(devTime,SID,sep="_"),SID)) %>% mutate(SID=recode(SID,"E3"="L8C","E4"="L16C"))
C14C19MC.ov <- read.delim("~/My_project/Extra_smncRNA/big_doc/sc_smallRNA_annotation/Human/C14MC_C19MC_cluster.ov.small.anno.bed",head=F,stringsAsFactors = F) %>% filter(V12=="miRNA") %>% select(V11,V4) %>% tbl_df() %>% unique() %>% rename(ID=V11,CMC=V4) 


temp.plot$human <- HE.counts.filter$mature %>% inner_join( C14C19MC.ov ,by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% right_join(temp.M %>% select(cell,SID),by="cell") %>% select(CMC,SID,cell,ID) %>% unique()%>% group_by(SID,cell,CMC) %>% summarise(nExpMiRNA=n_distinct(ID)) %>% right_join(C14C19MC.ov  %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")) %>%  mutate(SID=factor(SID,c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"),ordered = T)) %>% ggplot()+geom_boxplot(mapping=aes(x=CMC,y=nExpMiRNA/nTotalMiRNA,fill=SID),position="dodge",outlier.size = 0.25)+xlab("")+scale_fill_manual(values=c(lineage.col.set[c("L8CM","Prelineage","Early_ICM","Late_ICM","EB_TE","LB_TE")] %>% setNames(c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"))))+theme_classic()

cowplot::plot_grid(plotlist = temp.plot)
plot.results$CMC.perc <- temp.plot


if (check) {
  temp.M <-  ME.meta.filter%>% filter(batch!="Split" & EML!="unknown") %>% mutate(SID=EML) %>% mutate(SID=recode(SID,"L2and4C"="L2and4Cand8CandMorula","L8CM"="L2and4Cand8CandMorula"))
  C2C12MC.ov <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/C2MC_C12MC_ov_small_anno.bed",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V4 %in% c("C2MC","C12MC")) %>% filter(V12=="miRNA") %>% select(V4,V11,V12) %>% unique()
  
  ME.counts.filter$mature %>% inner_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11),by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% right_join(temp.M %>% select(cell,SID),by="cell")  %>% select(CMC,SID,cell,ID) %>% unique()%>% group_by(SID,cell,CMC) %>% summarise(nExpMiRNA=n_distinct(ID))  %>% right_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11) %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% group_by(SID,CMC,nTotalMiRNA) %>% summarise(nExpMiRNA=mean(nExpMiRNA)) 
}



if (FALSE) {
  temp.plot <- list()
  #' check the CMC one mouse
  temp.M <-  ME.meta.filter%>% filter(batch!="Split" & EML!="unknown") %>% mutate(SID=devTime)%>% mutate(devTime=recode(devTime,"64C"="MB","32C"="EB"))%>% mutate(SID=ifelse(EML %in% c("ICM","TE"),sub_EML,SID)) %>% mutate(SID=recode(SID,"Early_ICM"="EB_ICM","Late_ICM"="MB_ICM", "8C"="L8C","16C"="L16C","32C"="L32C","EB_TE"="Emergent_TE","MB_TE"="Early_TE","LB_TE"="Middle_TE"))  %>% filter(SID!="Emergent_TE") %>% mutate(SID=recode(SID,"Middle_TE"="MB_TE","Early_TE"="EB_TE"))
  C2C12MC.ov <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/C2MC_C12MC_ov_small_anno.bed",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V4 %in% c("C2MC","C12MC")) %>% filter(V12=="miRNA") %>% select(V4,V11,V12) %>% unique()
  

  temp.plot$mouse <- ME.counts.filter$mature %>% inner_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11),by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% right_join(temp.M %>% select(cell,SID),by="cell") %>% select(CMC,SID,ID) %>% unique()%>% group_by(SID,CMC) %>% summarise(nExpMiRNA=n_distinct(ID)) %>% right_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11) %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")) %>%  mutate(SID=factor(SID,c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"),ordered = T)) %>% ggplot()+geom_bar(mapping=aes(x=SID,y=nExpMiRNA/nTotalMiRNA,fill=CMC),stat="identity",position="dodge")+xlab("")
  
  
  
  
  #' check the CMC one human
  temp.M <-  HE.meta.filter %>% filter(batch!="HumanSplit" & devTime!="E7") %>% mutate(SID=EML) %>% mutate(devTime=recode(devTime,"E6"="MB","E7"="LB","E5"="EB"))%>% mutate(SID=ifelse(SID %in% c("ICM","TE"),paste(devTime,SID,sep="_"),SID)) %>% mutate(SID=recode(SID,"E3"="L8C","E4"="L16C"))
  C14C19MC.ov <- read.delim("~/My_project/Extra_smncRNA/big_doc/sc_smallRNA_annotation/Human/C14MC_C19MC_cluster.ov.small.anno.bed",head=F,stringsAsFactors = F) %>% filter(V12=="miRNA") %>% select(V11,V4) %>% tbl_df() %>% unique() %>% rename(ID=V11,CMC=V4) 


  temp.plot$human <- HE.counts.filter$mature %>% inner_join( C14C19MC.ov ,by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% right_join(temp.M %>% select(cell,SID),by="cell") %>% select(CMC,SID,ID) %>% unique()%>% group_by(SID,CMC) %>% summarise(nExpMiRNA=n_distinct(ID)) %>% right_join(C14C19MC.ov  %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")) %>%  mutate(SID=factor(SID,c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"),ordered = T)) %>% ggplot()+geom_bar(mapping=aes(x=SID,y=nExpMiRNA/nTotalMiRNA,fill=CMC),stat="identity",position="dodge")+xlab("")
  
  cowplot::plot_grid(plotlist = temp.plot)
  
  
  
  
  
}




#ME.counts.filter$mature %>% inner_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11),by="ID") %>% gather(cell,ct,-c(ID,CMC)) %>% filter(ct >0) %>% mutate(ct=1) %>% right_join(temp.M %>% select(cell,SID),by="cell") %>% group_by(ID,CMC,SID) %>% summarise(nExpCell=n_distinct(cell)) %>% inner_join(temp.M %>% group_by(SID) %>% summarise(nTotalCell=n_distinct(cell))) %>% group_by(SID,CMC) %>% summarise(nExpMiRNA=n_distinct(ID)) %>% right_join(C2C12MC.ov %>% rename(CMC=V4,ID=V11) %>% group_by(CMC) %>% summarise(nTotalMiRNA=n_distinct(ID)),by="CMC") %>% ggplot()+geom_point(mapping=aes(x=SID,y=nExpMiRNA/nTotalMiRNA,fill=CMC))

