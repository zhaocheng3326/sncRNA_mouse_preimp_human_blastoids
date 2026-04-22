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

#' average expression
miRNA.ave.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.ave.exp.rds"))


# taret paris
target.pairs <- readRDS(paste0("tmp_data/",TD,"/miRNA.gene.target.rds"))

#‘ check the seed freq 
print(
  miRNA.seed %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 1) %>% ggplot()+geom_bar(mapping = aes(x=reorder(seed,-nmiRNA),y=nmiRNA),stat="identity")+ theme_classic() + theme(axis.text.x=element_text(angle = 90))
)

miRNA.seed %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 1) %>% nrow()
miRNA.seed %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 2) %>% nrow()
miRNA.seed %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 3) %>% nrow()


#' only expressed
print(
  miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp)) %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 1) %>% ggplot()+geom_bar(mapping = aes(x=reorder(seed,-nmiRNA),y=nmiRNA),stat="identity")+ theme_classic() + theme(axis.text.x=element_text(angle = 90))
)
miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp))%>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 1) %>% nrow()
miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp))%>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 1)  %>% pull(nmiRNA) %>% sum()

miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp))%>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 2) %>% nrow()
miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp))%>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 3) %>% nrow()


#' check motif frequence in different type of miRNA
temp.plot <- list()
for (l in unique(miRNA.type$type)) {
 ln <-  miRNA.type %>% filter(type==l) %>% pull(miRNA) %>% unique() %>% length()
  temp.plot[[l]] <- miRNA.type %>% filter(type %in% l) %>% inner_join( miRNA.seed ,by="miRNA") %>% pull(seed) %>% ggseqlogo(seq_type='rna' )+ggtitle(paste0(l,"(",ln,")"))+FunTitle()+ylim(0,2)
} 
cowplot::plot_grid(plotlist = temp.plot,ncol=4)

#' check motif freq in different type of miRNA (only considering miRNA with seed freq >= 2)
temp.miRNA.seed <- miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp)) %>% group_by(seed) %>% mutate(nmiRNA=n_distinct(miRNA)) %>% filter(nmiRNA > 1) %>% unique()
temp.plot <- list()
for (l in unique(miRNA.type$type)) {
  ln <-  miRNA.type %>% filter(type==l) %>% pull(miRNA) %>% intersect(temp.miRNA.seed$miRNA) %>% unique() %>% length()
  if (ln > 4) {
    temp.plot[[l]] <- miRNA.type %>% filter(type %in% l) %>% inner_join(temp.miRNA.seed ,by="miRNA") %>% pull(seed) %>% ggseqlogo(seq_type='rna' )+ggtitle(paste0(l,"(",ln,")"))+FunTitle()+ylim(0,2)
  }
} 
cowplot::plot_grid(plotlist = temp.plot,ncol=4)

# check the purple one 
l="purple"
miRNA.type %>% filter(type %in% l) %>% inner_join( miRNA.seed ,by="miRNA")

# check the freq of seed distribution
freq.seed <- miRNA.seed %>% filter(miRNA %in% rownames(coseq.miRNA.exp)) %>% group_by(seed) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% arrange(desc(nmiRNA)) %>% filter(nmiRNA > 2) %>% pull(seed)


#' check the seed distribution 
#' lineage markers(seed)
miRNA.seed  %>% filter(seed %in% freq.seed ) %>% inner_join(miRNA.type %>% filter(!type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% group_by(seed,type) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% spread(type,nmiRNA) %>% as.data.frame() %>% tibble::column_to_rownames("seed") %>% replace(.,is.na(.),0) %>% pheatmap(scale="none",display_numbers = T,number_format="%.0f",main="Number of miRNAs")


#' modules (seed)
miRNA.seed  %>% filter(seed %in% freq.seed ) %>% inner_join(miRNA.family %>% rename(miRNA=mature_miRNA),by="miRNA") %>% mutate(SID=paste(seed)) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% group_by(SID,type) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% spread(type,nmiRNA) %>% as.data.frame() %>% tibble::column_to_rownames("SID") %>% replace(.,is.na(.),0) %>% pheatmap(scale="none",display_numbers = T,number_format="%.0f",main="Number of miRNAs")

#' check AAGUGCU
temp.sel.seed <- "AAGUGCU"
temp <- miRNA.seed  %>% filter(seed %in% temp.sel.seed  ) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 

temp %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,type) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,annotation_colors = list(type=c("brown"="brown","turquoise"="turquoise")),main=temp.sel.seed)


#' check 
temp.sel.seed <- c("AGCAGCA","AUUGCAC","ACCCGUA","GUAAACA","AAAGUGC","AGUGCAA")
temp <- miRNA.seed  %>% filter(seed %in% temp.sel.seed  ) %>% inner_join(miRNA.type %>% filter(!type %in% miRNA.MEs.anno$miRNA_mdc ) %>% select(miRNA),by="miRNA")  %>% left_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 

temp %>% arrange(seed) %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,seed) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,cluster_rows=T,main=paste(temp.sel.seed,collapse=","))


#' check ACAUUCA & GAGGUAG in  oocyte markers
temp.sel.seed <- c("ACAUUCA","GAGGUAG")
temp <- miRNA.seed  %>% filter(seed %in% temp.sel.seed  ) %>% inner_join(miRNA.type %>% filter(type %in% c("oocyte")) %>% select(miRNA),by="miRNA")  %>% left_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 

temp %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,seed) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,main=paste(temp.sel.seed,collapse=","))


#' modules (seed+family)
miRNA.seed  %>% filter(seed %in% freq.seed ) %>% inner_join(miRNA.family %>% rename(miRNA=mature_miRNA),by="miRNA") %>% mutate(SID=paste(seed,miRNA_family,sep="_")) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% group_by(SID,type) %>% summarise(nmiRNA=n_distinct(miRNA)) %>% spread(type,nmiRNA) %>% as.data.frame() %>% tibble::column_to_rownames("SID") %>% replace(.,is.na(.),0) %>% pheatmap(scale="none",display_numbers = T,number_format="%.0f",main="NO.of miRNAs")

#' check
temp.sel.sid <- c("GAGGUAG_let-7:MIPF0000002")
temp <- miRNA.seed  %>% filter(seed %in% freq.seed ) %>% inner_join(miRNA.family %>% rename(miRNA=mature_miRNA),by="miRNA") %>% mutate(SID=paste(seed,miRNA_family,sep="_")) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% filter(SID %in% temp.sel.sid) %>% arrange(type) %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 

temp %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,type) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,annotation_colors = list(type=c("green"="green","magenta"="magenta","purple"="purple")),main=temp.sel.sid)

target.pairs %>% filter(miRNA %in% c("mmu-let-7e-5p","mmu-let-7c-5p","mmu-let-7a-5p","mmu-let-7d-5p")) %>% split(.,.$miRNA) %>% lapply(function(x){x$gene}) %>% ggvenn::ggvenn()#,"mmu-let-7f-5p","mmu-let-7g-5p","mmu-let-7i-5p","mmu-miR-98-5p"

data.temp <- subset(readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds")),cell=(meta.filter %>% filter(batch!="Split" & EML!="unknown") %>% pull(cell)))
data.temp@meta.data$EML <- factor((meta.filter %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"],EML.od,ordered = T)
VlnPlot(data.temp,temp$miRNA,group.by = "EML")


#' check
temp.sel.sid <- c("GUGUGCA_mir-467:MIPF0000316")
temp <- miRNA.seed  %>% filter(seed %in% freq.seed ) %>% inner_join(miRNA.family %>% rename(miRNA=mature_miRNA),by="miRNA") %>% mutate(SID=paste(seed,miRNA_family,sep="_")) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% filter(SID %in% temp.sel.sid) %>% arrange(type) %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 
temp %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,type) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,annotation_colors = list(type=c("blue"="blue","turquoise"="turquoise")),main=temp.sel.sid)
target.pairs %>% filter(miRNA %in% c("mmu-let-7e-5p","mmu-let-7c-5p","mmu-let-7a-5p","mmu-let-7d-5p")) %>% mutate(miRNA=gsub("mmu-","",miRNA)) %>% split(.,.$miRNA) %>% lapply(function(x){x$gene}) %>% ggvenn::ggvenn(text_size=3,set_name_size=4)#,"mmu-let-7f-5p","mmu-let-7g-5p","mmu-let-7i-5p","mmu-miR-98-5p"



#' check AAGUGCU
temp.sel.seed <- "AAGUGCU"
temp <- miRNA.seed  %>% filter(seed %in% temp.sel.seed  ) %>% inner_join(miRNA.type %>% filter(type %in% miRNA.MEs.anno$miRNA_mdc ),by="miRNA") %>% inner_join(miRNA.ave.exp %>% rename(miRNA=gene),by="miRNA") 
temp %>% select(miRNA,all_of(EML.od)) %>% tibble::column_to_rownames("miRNA") %>% pheatmap(scale="row",annotation_row = (temp %>% select(miRNA,type) %>% tibble::column_to_rownames("miRNA")),cluster_cols = F,annotation_colors = list(type=c("brown"="brown","turquoise"="turquoise")),main=temp.sel.seed)
miRNA.family %>% filter(mature_miRNA %in% temp$miRNA)
