#' ---
#' title: answer mimic abundance
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

###The Dlk1-Dio3 region encodes an unusually large number of miRNA species, so counting species alone can artificially inflate apparent enrichment relative to other loci. Without a null model (e.g., permutation preserving family structure), the 79.5% value reported in Figure 3a,c may reflect cluster size rather than true ICM specificity. Support this claim beyond descriptive visualization requires a statistical enrichment test and abundance-containing metrics.

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


# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=ifelse(EML=="prelineage" & devTime %in% c("2C","4C"),"L2and4C",EML)) %>% mutate(EML=ifelse(EML=="prelineage" & !devTime %in% c("2C","4C"),"L8CM",EML))  %>% select(cell,RNA_EML,small_EML,EML) #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
batch4.meta.filter <- meta.filter %>% filter(batch=="batch4" ) 

meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

sel.samples <- c("Sample_1558_0","Sample_1557_0","Sample_1556_0","Sample_1555_0","Sample_1554_0","Sample_1553_0","Sample_1619_0","Sample_1618_0","Sample_1617_0","Sample_1616_0","Sample_1615_0","Sample_1614_0","Sample_1613_0","Sample_1612_0")

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))

temp.plot <- list()

temp.raw <- batch4.meta.filter %>% filter(cell %in% sel.samples) %>% mutate(SID=recode(stage,"10 uM mimic ctrl 8 cell"="Control","10 uM novel miR 8 cell"="Treated")) %>% select(cell,SID) %>% left_join(read.delim("tmp_data/batch4/check_mimic/mimic.step1.umi.stat",head=F) %>% tbl_df() %>% setNames(c("cell","UMI")),by="cell") %>% replace(.,is.na(.),0) %>% mutate(miRNA="miR-29031-3p")
temp.exp <- temp.raw%>% group_by(miRNA,SID) %>% summarise(mean=mean(UMI),sd=sd(UMI),n=n_distinct(cell)) %>% mutate(sem=sd/(n^0.5)) %>% mutate(SID=factor(SID,c("Control","Treated"),ordered =T))

g="miR-29031-3p"
temp.plot[[paste0(g,"_umi")]] <- temp.exp %>% filter(miRNA %in% g)%>%  ggplot(mapping=aes(x=SID,y=mean,fill=SID))+geom_bar(stat="identity",position="dodge",width=0.75)+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),size=0.25,width=0.5,position=position_dodge(1))+geom_point(data=temp.raw%>% filter(miRNA %in% g),mapping=aes(x=SID,y=UMI),position = position_jitter(width = 0.2,height = 0))+theme_classic()+xlab("")+ylab("NO. of UMI")+ggtitle(g)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=c("Control"="grey88","Treated"="grey33"))+NoLegend()


#mmu-miR-221-3p
temp.raw <- data.all.ob.umap %>% filter(sub_EML %in% c("EB_TE")) %>% mutate(SID="Emergent_TE") %>% bind_rows(data.all.ob.umap %>% filter(EML=="L8CM" & devTime=="16C") %>% mutate(SID="Morula")) %>% select(cell,SID) %>% inner_join(counts.filter$mature %>% filter(ID=="mmu-miR-221-3p") %>% gather(cell,UMI,-ID),b="cell") %>% rename(miRNA=ID) 
temp.exp <- temp.raw%>% group_by(miRNA,SID) %>% summarise(mean=mean(UMI),sd=sd(UMI),n=n_distinct(cell)) %>% mutate(sem=sd/(n^0.5)) %>% mutate(SID=factor(SID,c("Morula","Emergent_TE"),ordered =T))

g="mmu-miR-221-3p"
temp.plot[[paste0(g,"_EE")]] <- temp.exp %>% filter(miRNA %in% g)%>%  ggplot(mapping=aes(x=SID,y=mean,fill=SID))+geom_bar(stat="identity",position="dodge",width=0.75)+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),size=0.25,width=0.5,position=position_dodge(1))+theme_classic()+xlab("")+ylab("NO. of UMI")+ggtitle(g)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=c("Morula"="#C49A00","Emergent_TE"="#B3E5E1"))+NoLegend()+ylim(0,7.5)#+geom_point(data=temp.raw%>% filter(miRNA %in% g) %>% mutate(UMI=ifelse(UMI >20,20,UMI)),mapping=aes(x=SID,y=UMI),position = position_jitter(width = 0.2,height = 0))
cowplot::plot_grid(plotlist = temp.plot)

load(paste0("tmp_data/",TD,"/lineage.segregation.DEG.out.Rdata"),verbose=T)
DEG.results.list$miRNA$EB_TE_vs_L16C$DEG.result %>% filter(gene==g)
pdf("tmp_data/temp.R1Q4.pdf",4.5,4)
cowplot::plot_grid(plotlist = temp.plot)
dev.off()
