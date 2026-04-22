#' ---
#' title: check the enrichment of seed sequence
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
  library(ggseqlogo)
  library(ComplexUpset)
  #library(ff)
  #library(scran)
  #library(WGCNA)
  #library(hdWGCNA)
  #library(UCell)
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

#' miRNA family information
miRNA.family <- readRDS(paste0("tmp_data/",TD,"/miRNA.family.rds"))

#' seed sequence information
miRNA.seed <- read.delim("big_doc/sc_smallRNA_annotation/Mouse/mmu.mature.miRNA.seed.seq.tsv",stringsAsFactors = F,head=F) %>% tbl_df() %>% setNames(c("miRNA","seed")) %>% unique()

#' check the seed distribution within each miRNA family
temp <- miRNA.family %>% inner_join(miRNA.seed %>% rename(mature_miRNA=miRNA) ,by="mature_miRNA") %>% unique() 
temp %>% group_by(miRNA_family,seed) %>% summarise(nM=n_distinct(mature_miRNA))
temp %>% filter(grepl("-5p",mature_miRNA)) %>% group_by(miRNA_family) %>% summarise(nS=n_distinct(seed)) %>% arrange(desc(nS))
temp %>% filter(grepl("-5p",mature_miRNA)) %>% filter(miRNA_family=="mir-467:MIPF0000316") %>% pull(seed) %>% unique() %>% length()



#' check the miRNA family distribution shared each seed
temp <- miRNA.family %>% inner_join(miRNA.seed %>% rename(mature_miRNA=miRNA) ,by="mature_miRNA") %>% unique() 
temp %>% group_by(seed) %>% mutate(nMF=n_distinct(miRNA_family),nM=n_distinct(mature_miRNA)) %>% filter(nMF >1) %>% select(-mature_miRNA) %>% unique() %>% arrange(desc(nMF),seed)


temp %>% filter(grepl("-5p",mature_miRNA)) %>% group_by(miRNA_family) %>% summarise(nS=n_distinct(seed)) %>% arrange(desc(nS))
temp %>% filter(grepl("-5p",mature_miRNA)) %>% filter(miRNA_family=="mir-467:MIPF0000316")
temp %>% group_by(seed,miRNA_family) %>% mutate(nM=n_distinct(mature_miRNA)) %>% filter(nM >1) %>% ungroup() %>%  select(miRNA_family,seed) %>% unique()




#' loading information
TF_miRNA <- readRDS(paste0("tmp_data/",TD,"/TF_miRNA.rds"))
miRNA.bed <- readRDS(paste0("tmp_data/",TD,"/miRNA.bed.rds"))
miRNA.mk.out <- readRDS(file=paste0("tmp_data/",TD,"/","allCells.miRNA.mk.rds"))


# all cells  annotation
coseq.small.umap <- readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% mutate(EML=recode(EML,"prelineage"="L8CM")) %>% select(cell,RNA_EML,small_EML,EML)  #%>% select(cell,EML,RNA_EML,small_EML)
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% mutate(sub_EML=EML) %>% rows_update(coseq.small.umap %>% filter(EML!="unknown") %>% select(cell,EML) %>% mutate(sub_EML="None"),by="cell")  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch)
meta.filter <- meta.filter %>% filter(! batch %in% c("batch4","batch5")) %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))
#' update the full cell annotation
meta.filter <- meta.filter  %>% inner_join(data.all.ob.umap %>% select(cell,sub_EML,EML),by="cell")
#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(RNA_EML="None",small_EML="None") %>% rows_update(coseq.small.umap %>% filter(EML!="unknown"),by="cell")%>% mutate(batch=recode(batch,"Split1"="Split","Split2"="Split"))


#' coseq smt2 part
coseq.gene.exp <- readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.gene.exp.log.rds"))
coseq.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds")) ### including other batches


#' loading ME for coseq
md.miRNA.list <- readRDS(paste0("tmp_data/",TD,"/TOM/","coseq_proj",".data.query.perum.md.list.rds"))
miRNA.MEs <-md.miRNA.list$MEs[,colnames(md.miRNA.list$MEs) %>% setdiff(c("grey"))]
miRNA.MEs.anno <- md.miRNA.list$md.anno
colnames(miRNA.MEs ) <- (miRNA.MEs.anno %>% tibble::column_to_rownames("miRNA_mdc"))[colnames(miRNA.MEs ),"miRNA_md"]

#' using the cells have the same RNA_EML and small_EML
sel_cells <- meta.filter %>% filter(RNA_EML==small_EML & batch=="Split") %>% pull(cell) %>% intersect(readRDS(paste0("tmp_data/",TD,"/mouse_coseq_smt2.updated.meta.filter.rds"))$cell)%>% intersect(rownames(md.miRNA.list$MEs))
coseq.miRNA.exp <- coseq.miRNA.exp[,sel_cells]

#' loading traj related miRNAs
load(paste0("tmp_data/",TD,"/mouse.miRNA.SCP.main.traj.Rdata"),verbose = T)

#' miRNA types
miRNA.type <- md.miRNA.list$modules %>% rename(miRNA=gene_name,type=module) %>% select(miRNA,type) %>% tbl_df() %>% filter(type!="grey") %>%  bind_rows(miRNA.mk.out %>% rename(miRNA=gene,type=cluster) %>% select(miRNA,type) ) %>% bind_rows(psdt.genes$ICM_main_traj %>% mutate(type=paste0("ICM_traj_",cluster),miRNA=gene) %>% select(miRNA,type)) %>% bind_rows(psdt.genes$TE_main_traj %>% mutate(type=paste0("TE_traj_",cluster),miRNA=gene) %>% select(miRNA,type))

#' normalized expression
miRNA.norm.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

#' average expression
miRNA.ave.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.ave.exp.rds"))

# taret pairs
target.pairs <- readRDS(paste0("tmp_data/",TD,"/miRNA.gene.target.rds"))

#' counts filter 
counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))


#' loading gene meta data
load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)



cell.list <- meta.filter %>% filter(sub_EML %in% c("Early_ICM","EB_TE","Late_ICM","MB_TE","LB_TE")) %>% filter(batch %in% c("batch1","batch2","batch4")) %>% mutate(SID=sub_EML) %>% bind_rows(meta.filter %>% filter(EML!="unknown") %>% filter(batch %in% c("batch1","batch2","batch4")) %>% mutate(SID=EML)) %>% select(cell,SID) %>% unique() %>% split(.,.$SID) %>% lapply(function(x){x$cell})
cell.list$ICM_TE <- c(cell.list$ICM,cell.list$TE)
cell.list %>% lapply(length) %>% unlist() %>% min()

#' check the rank
#' check the abundance of mmu-miR-221-3p in EB_TE cells
temp.cells <- meta.filter %>% filter(batch!="Split") %>% filter(sub_EML=="EB_TE") %>% pull(cell)
temp.exp <- miRNA.norm.exp[,temp.cells] %>% as.data.frame()

temp.rank <- temp.exp %>% tibble::rownames_to_column("ID") %>% gather(cell,logExp,-ID) %>% tbl_df() %>% arrange(cell,desc(logExp)) %>% split(.,.$cell) %>% lapply(function(x){x %>%tibble::rowid_to_column("rank")}) %>% do.call("bind_rows",.) %>% filter(ID=="mmu-miR-221-3p") %>% pull(rank) 
summary(temp.rank)
temp.valid.rank <- median(temp.rank)
print(temp.valid.rank)

#" check the proportion of special miRNA in different cell groups
temp.out <- list()
for (g in names(cell.list)) {
  temp <-  miRNA.norm.exp[,cell.list[[g]]] %>% as.data.frame()%>% tibble::rownames_to_column("ID") %>% tbl_df() %>% gather(cell,logExp,-ID) %>% tbl_df() %>% arrange(cell,desc(logExp)) %>% split(.,.$cell) %>% lapply(function(x){x %>%tibble::rowid_to_column("rank")}) %>% do.call("bind_rows",.) 
  temp.out[[g]] <- temp %>% group_by(ID) %>% summarise(rank.median=median(rank)) %>% mutate(SID=g)
}
temp.out <- temp.out %>% do.call("bind_rows",.)
temp.out %>% rename(miRNA=ID) %>% inner_join(miRNA.type,by="miRNA",relationship = "many-to-many") %>% filter(rank.median < temp.valid.rank) %>% group_by(SID,type) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% inner_join(miRNA.type %>% group_by(type) %>% summarise(total_miRNA=n_distinct(miRNA)),by="type")  %>% mutate(prop=nmiRNA/total_miRNA) %>% select(SID,type,prop) %>% spread(type,prop)%>% replace(.,is.na(.),0) %>% tibble::column_to_rownames("SID") %>% pheatmap(display_numbers = T,main="Proportion of cells passed the cutoff(miR-221-3p)")


#" check the proportion of special miRNA in different cell groups
temp.plot <- list()
for (g in names(cell.list)) {
  temp <-  miRNA.norm.exp[,cell.list[[g]]] %>% as.data.frame()%>% tibble::rownames_to_column("ID") %>% tbl_df() %>% gather(cell,logExp,-ID) %>% tbl_df() %>% arrange(cell,desc(logExp)) %>% split(.,.$cell) %>% lapply(function(x){x %>%tibble::rowid_to_column("rank")}) %>% do.call("bind_rows",.)  %>% rename(miRNA=ID) %>% inner_join(miRNA.type,by="miRNA",relationship = "many-to-many")
  for (m in c(c(EML.od,"ICM_TE")) ) {
    if(m==g) {
      temp.plot[[paste0("cg:",g," mg:",m)]] <- temp %>% filter(type==m) %>% ggplot()+geom_histogram(mapping=aes(x=rank),fill="royalblue3",bins=20)+theme(plot.title = element_text(hjust=0.5))+geom_vline(xintercept = temp.valid.rank,linetype="dashed")+theme_classic()+ylab("No. of cells")+ggtitle(paste0("cg:",g," mg:",m))+FunTitle()
    }
  }
}
cowplot::plot_grid(plotlist = temp.plot)
  
