#' ---
#' title: leave one embryo l1 model as reviewer suggested
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
  #library(scran)
  #library(batchelor)
  #library(SeuratWrappers)
  #library(slingshot)
  #library(tradeSeq)
  library(glmnet)
  #library(pROC)
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


#' loading gene meta data
load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))

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

#' select need cells
meta.filter <- meta.filter %>% filter(batch %in% c("batch1","batch2","batch4"))  %>% filter(EML %in% c("ICM","TE"))


#' loading expression
sel.exp.list <- list()
sel.exp.list<- readRDS(paste0("tmp_data/",TD,"/","allCells.other.sncRNA",".norm.list.rds"))
sel.exp.list$miRNA <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))



#' select test embryo with at least 5 cells, contain ICM & TE, otherwise unable to calculate the auc values
#meta.filter %>% group_by(embryo,EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(EML,nCell)
#sel.embryo <- meta.filter %>% group_by(embryo) %>% summarise(nCell=n_distinct(cell)) %>% filter(nCell>=5) %>% pull(embryo)
sel.embryo <- meta.filter %>% group_by(embryo,EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(EML,nCell) %>% filter(!is.na(ICM) & !is.na(TE)) %>% filter((ICM+TE) >=5) %>% pull(embryo) ## no obivious differene setting to 20


#setting
#sel.type <- c("miRNA","snoRNA")
sel.type <- c("miRNA","snoRNA","tRNA","snRNA","piRNA","rRNA")
zs.limit <- 2
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

plot.results <- list()

if (file.exists(paste0("tmp_data/",TD,"/","L1.regression.Rdata"))) {
  load(paste0("tmp_data/",TD,"/","L1.regression.Rdata"),verbose=T)
}else{
  temp.perf.list <- list()
  temp.coef.list <- list()
  for (st in sel.type) {
    tag <- paste("real",st,sep=":")
    temp.out <- Fun_ICM_TE_L1_regression(sel.exp.list[[st]],meta.filter,sel.embryo,tag) 
    temp.perf.list[[st]] <- temp.out$perf
    temp.coef.list[[st]] <- temp.out$coef
  }

  perf.out <- temp.perf.list %>% do.call("bind_rows",.) %>% tbl_df()
  coef.out <- temp.coef.list %>% do.call("bind_rows",.) %>% tbl_df()
  save(perf.out,coef.out,file=paste0("tmp_data/",TD,"/","L1.regression.Rdata"))
}

plot.results$L1.ICM_TE.perf <- perf.out %>% mutate(type=gsub("real:","",tag)) %>% select(-tag) %>% gather(feature,score,-c(sel_embryo,type)) %>% mutate(type=factor(type,sel.type,ordered = T)) %>% filter(feature %in% c("auc","ba_acc")) %>% mutate(feature=recode(feature,"auc"="AUC","ba_acc"="balanced accuracy"))  %>% ggplot(mapping=aes(x=feature,y=score))+geom_boxplot(mapping=aes(fill=type),color="grey50",width=.5,coef = 1.5,outlier.size=0.2)+ theme_classic()  +theme(axis.text.x=element_text(angle = 30,hjust = 1),panel.background = element_rect(fill = NA,colour="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.grid.major.x = element_line(colour = NA),strip.text = element_text(face ="bold"),strip.background = element_blank()) +scale_fill_manual(values=type.col.set)+ggtitle("L1 (ICM vs TE)")+FunTitle()+xlab("")+ylab("Value")
plot.results$L1.ICM_TE.perf 



#’ check top ranking features
temp.top.rank.features <- coef.out %>% mutate(type=gsub("real:","",tag)) %>% filter(type!="snRNA")%>% select(-tag) %>% group_by(ID,type) %>% summarise(coef=median(coef)) %>% mutate(UpDown=ifelse(coef< 0,"ICM","TE")) %>% group_by(type,UpDown) %>% top_n(10,abs(coef))

temp.plot <- list()
for (st in sel.type) {
  if (st!="snRNA") {
    temp.plot[[st]] <- temp.top.rank.features %>% filter(type==st) %>% ggplot()+geom_bar(mapping=aes(y=reorder(ID,coef),x=coef),stat="identity")+theme_classic()+ggtitle(st)+FunTitle()+ylab("")+xlab("coefficient")
  }
}
cowplot::plot_grid(plotlist = temp.plot)

#' check the expressed percentage of those top ranks
temp.top.rank.features <- coef.out %>% mutate(type=gsub("real:","",tag)) %>% filter(type!="snRNA")%>% select(-tag) %>% group_by(ID,type) %>% summarise(coef=median(coef)) %>% mutate(UpDown=ifelse(coef< 0,"ICM","TE")) %>% group_by(type,UpDown) %>% top_n(10,abs(coef))

print(
  counts.filter$mature %>% inner_join(temp.top.rank.features %>% ungroup()%>% select(-UpDown),by="ID") %>% gather(cell,ct,-c(ID,type,coef)) %>% inner_join(meta.filter %>% select(cell,EML),by="cell") %>% filter(ct > 0) %>% group_by(ID,type,EML) %>% summarise(nCell=n_distinct(cell)) %>% inner_join(meta.filter %>% select(cell,EML) %>% group_by(EML) %>% summarise(nBGcell=n_distinct(cell)),by="EML") %>% mutate(prop=nCell/nBGcell) %>% mutate(type=factor(type,sel.type,ordered = T)) %>% ggplot()+geom_boxplot(mapping=aes(x=type,fill=EML,y=prop))+ylab("Proportion of expressed cells")+xlab("")+ylim(0,1)+theme_classic()
)
plot.results$miRNA.top.rf <- counts.filter$mature %>% inner_join(temp.top.rank.features %>% ungroup()%>% select(-UpDown),by="ID") %>% gather(cell,ct,-c(ID,type,coef)) %>% inner_join(meta.filter %>% select(cell,EML),by="cell") %>% filter(ct > 0) %>% group_by(ID,type,EML) %>% summarise(nCell=n_distinct(cell)) %>% inner_join(meta.filter %>% select(cell,EML) %>% group_by(EML) %>% summarise(nBGcell=n_distinct(cell)),by="EML") %>% mutate(prop=nCell/nBGcell) %>% mutate(type=factor(type,sel.type,ordered = T)) %>% filter(type=="miRNA") %>% ggplot()+geom_boxplot(mapping=aes(x=EML,fill=EML,y=prop),width=0.6)+ylab("Proportion of expressed cells")+xlab("")+theme_classic()+ylim(0,0.05)
plot.results$miRNA.top.rf



#‘ sample permutation
#' generate sample meta.filter
if (file.exists(paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.meta.rds"))) {
  perm.meta.filter.list <- readRDS(paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.meta.rds"))
}else{
  # must save to repeat
  perm.meta.filter.list <- list()
  for (n in 1:1000) {
    temp.meta.filter <- meta.filter %>% select(cell,embryo,EML) %>% mutate(perm=paste0("perm_",n))
    temp.meta.filter$EML=sample(temp.meta.filter$EML)
    perm.meta.filter.list[[paste0("perm_",n)]] <- temp.meta.filter
  }
  saveRDS(perm.meta.filter.list,paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.meta.rds"))
}

# permutation calculation
if (file.exists(paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.perm.perf.rds"))) {
  perm.perf.out <- readRDS(paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.perm.perf.rds"))
}else{
  perm.perf.list <- list()
  #perm.coef.list <- list()

  for (st in sel.type) {
    temp.list <- list()
    for (n in 1:1000) {
      temp.list[[paste0("perm_",n)]] <- paste0("perm_",n)
    }
    perm.perf.list[[st]] <- foreach (x=temp.list,n=names(temp.list),.combine=c) %dopar% {
      rv=list()
      rv[[n]]=Fun_para_ICM_TE_L1_regression(sel.exp.list,perm.meta.filter.list,x,sel.embryo,st)
      rv
    }
  }
  
  
  perm.perf.out <- perm.perf.list %>% lapply(function(x) {do.call("bind_rows",x)}) %>% do.call("bind_rows",.) %>% tbl_df()
  saveRDS(perm.perf.out,paste0("tmp_data/",TD,"/","L1.regression.perm.ICM_TE.perm.perf.rds"))
}
#<- meta.filter %>% select(cell,EML)
perm.perf.out %>% separate(tag,c("tag","type"),sep=":") %>% gather(feature,score,-c(sel_embryo,type,tag)) %>% group_by(type,feature,tag) %>% summarise(perm_score=median(score)) %>%ggplot()+ geom_histogram(mapping=aes(x=perm_score))+facet_grid(feature~type)
#' get pvalues 
perm.perf.out %>% separate(tag,c("tag","type"),sep=":") %>% gather(feature,score,-c(sel_embryo,type,tag)) %>% group_by(type,feature,tag) %>% summarise(perm_score=median(score)) %>% inner_join(perf.out %>% mutate(type=gsub("real:","",tag)) %>% select(-tag) %>% gather(feature,score,-c(sel_embryo,type)) %>% group_by(type,feature) %>% summarise(real_score=median(score)),by=c("type","feature")) %>% filter(perm_score < real_score) %>% group_by(type,feature) %>% summarise(pvalue=1-n_distinct(tag)/1000) %>% spread(feature,pvalue)




