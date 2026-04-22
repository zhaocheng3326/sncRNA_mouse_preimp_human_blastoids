#' ---
#' title: cross-species"
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




#' get expression data
if (file.exists(paste0("tmp_data/",TD,"/","species.miRNA.ave.exp.rds"))) {
  species.miRNA.exp <- readRDS(paste0("tmp_data/",TD,"/","species.miRNA.ave.exp.rds"))
}else{
  species.miRNA.exp <- list()
  
  #' loading sorted number of genes
  miRNA.sort.nExp.miRNA.list <- readRDS(file=paste0("tmp_data/",TD,"/","miRNA.sort.nExp.miRNA.list.rds"))

  #' human data
  load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T) 
  data.human.umap <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.coord.check.rds")
  human.norm.exp <- readRDS(paste0("~/My_project/Extra_smncRNA/tmp_data/","Apr_2023","/","miRNA",".norm.rds"))$mBN
  

  species.miRNA.exp$human.exp.miRNA <- miRNA.sort.nExp.miRNA.list$human.E3E6
  
  data.temp.umap <- data.human.umap %>% filter(batch!="HumanSplit" & devTime!="E7") %>% mutate(SID=recode(EML,"EarlyTE"="TE","mural"="TE","polar"="TE","E3"="L8CM","E4"="L8CM")) %>% filter(SID %in% c("ICM","TE","L8CM"))
  #data.temp.umap <- data.human.umap %>% filter(batch!="HumanSplit") %>% mutate(SID=recode(EML,"EarlyTE"="TE","mural"="TE","polar"="TE","E3"="L8CM","E4"="L8CM")) %>% filter(SID %in% c("ICM","TE","L8CM"))
  temp.counts <- (readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")$mature %>% tibble::column_to_rownames("ID"))[rownames(human.norm.exp),data.temp.umap$cell]
  data.sub.temp <- CreateSeuratObject(temp.counts, meta.data = (data.temp.umap %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.sub.temp@assays$RNA$data <- as.sparse(human.norm.exp[,colnames(data.sub.temp)])
  Idents(data.sub.temp) <- factor(data.sub.temp@meta.data$SID)
  species.miRNA.exp$human.t3 <- AverageExpression(data.sub.temp,group.by = "SID")$RNA %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% gather(cluster,ave_exp,-gene) %>% tbl_df() %>% mutate(species="human") %>% mutate(cluster=gsub("-","_",cluster)) %>% rename(SID=cluster)
  
  species.miRNA.exp$human.t3.mk <- FindAllMarkers(data.sub.temp) %>% tbl_df()%>% filter(p_val_adj < 0.05)  
  species.miRNA.exp$human.t3.mk.hl <- species.miRNA.exp$human.t3.mk %>% filter(avg_log2FC >0) %>% group_by(cluster) %>% top_n(5,-1*p_val_adj) %>% select(gene,cluster) %>% unique() %>% pull(gene)
  #species.miRNA.exp$human.t3.mk <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.FM.mk.rds")$bigL %>% filter(p_val_adj < 0.05) %>% filter(cluster!="unclassified") %>%select(gene,cluster) %>% unique()
  #species.miRNA.exp$human.t3.mk.hl <- readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/miRNA.IT.FM.mk.rds")$bigL %>% filter(p_val_adj < 0.01) %>% filter(cluster!="unclassified") %>% group_by(cluster) %>% top_n(7,-1*p_val_adj) %>% select(gene,cluster) %>% unique() %>% pull(gene)
  
  data.temp.umap <- data.human.umap %>% filter(batch!="HumanSplit" & devTime!="E7") %>% mutate(SID=recode(EML,"EarlyTE"="TE","mural"="TE","polar"="TE")) %>% mutate(devTime=recode(devTime,"E6"="MB","E7"="LB","E5"="EB"))%>% mutate(SID=ifelse(SID %in% c("ICM","TE"),paste(devTime,SID,sep="_"),SID)) %>% mutate(SID=recode(SID,"E3"="L8C","E4"="L16C"))
  temp.counts <- (readRDS("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/small.counts.filter.rds")$mature %>% tibble::column_to_rownames("ID"))[rownames(human.norm.exp),data.temp.umap$cell]
  data.sub.temp <- CreateSeuratObject(temp.counts, meta.data = (data.temp.umap %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.sub.temp@assays$RNA$data <- as.sparse(human.norm.exp[,colnames(data.sub.temp)])
  data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  species.miRNA.exp$human.dev <- AverageExpression(data.sub.temp,group.by = "SID")$RNA %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% gather(cluster,ave_exp,-gene) %>% tbl_df() %>% mutate(species="human")%>% mutate(cluster=gsub("-","_",cluster))%>% rename(SID=cluster)
  
  
  #' mouse
  data.mouse.umap <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% filter(!devTime %in% c("2C","4C","oocyte","sperm")) %>% mutate(sub_EML=EML) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".sub.data.ob.umap.rds")) %>% select(cell,sub_EML),by="cell")
  data.temp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.rds"))
  species.miRNA.exp$mouse.exp.miRNA <- miRNA.sort.nExp.miRNA.list$mouse.8C_64C
  
  data.temp.umap <- data.mouse.umap %>% filter(batch!="Split"& EML!="unknown") %>% mutate(SID=EML) %>% filter(SID %in% c("ICM","TE","L8CM"))
  data.sub.temp <- subset(data.temp,cells=data.temp.umap$cell)
  data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  Idents(data.sub.temp) <- factor(data.sub.temp@meta.data$SID)
  species.miRNA.exp$mouse.t3 <- AverageExpression(data.sub.temp,group.by = "SID")$RNA %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% gather(cluster,ave_exp,-gene) %>% tbl_df() %>% mutate(species="mouse") %>% mutate(cluster=gsub("-","_",cluster))%>% rename(SID=cluster)
  species.miRNA.exp$mouse.t3.mk <- FindAllMarkers(data.sub.temp) %>% tbl_df()%>% filter(p_val_adj < 0.05)  
  species.miRNA.exp$mouse.t3.mk.hl <- species.miRNA.exp$mouse.t3.mk %>% filter(avg_log2FC >0) %>% group_by(cluster) %>% top_n(5,-1*p_val_adj) %>% select(gene,cluster) %>% unique() %>% pull(gene)
  
  #miRNA.mk.out <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA.FM.mk.rds"))
  #species.miRNA.exp$mouse.t3.mk <- miRNA.mk.out$lineage %>% filter(p_val_adj < 0.05)  %>% filter(cluster %in% c("L8CM","TE","ICM")) %>% bind_rows(miRNA.mk.out$early_late %>% filter(p_val_adj < 0.05)) %>% select(gene,cluster) %>% unique()
  #species.miRNA.exp$mouse.t3.mk.hl <- miRNA.mk.out$lineage %>% filter(p_val_adj < 0.001)  %>% filter(cluster %in% c("L8CM","TE","ICM")) %>% bind_rows(miRNA.mk.out$early_late %>% filter(p_val_adj < 0.001)) %>% group_by(cluster) %>% top_n(7,-1*p_val_adj) %>% select(gene,cluster) %>% unique() %>% pull(gene)
  
  data.temp.umap <- data.mouse.umap %>% filter(batch!="Split" & EML!="unknown") %>% mutate(SID=devTime)%>% mutate(devTime=recode(devTime,"64C"="MB","32C"="EB"))%>% mutate(SID=ifelse(EML %in% c("ICM","TE"),sub_EML,SID)) %>% mutate(SID=recode(SID,"Early_ICM"="EB_ICM","Late_ICM"="MB_ICM", "8C"="L8C","16C"="L16C","32C"="L32C","EB_TE"="Emergent_TE","MB_TE"="Early_TE","LB_TE"="Middle_TE"))  %>% filter(SID!="Emergent_TE") %>% mutate(SID=recode(SID,"Middle_TE"="MB_TE","Early_TE"="EB_TE")) ## fix the name issues
  data.sub.temp <- subset(data.temp,cells=data.temp.umap$cell)
  data.sub.temp@meta.data$SID <- (data.temp.umap %>% tibble::column_to_rownames("cell"))[rownames(data.sub.temp@meta.data),"SID"]
  species.miRNA.exp$mouse.dev <- AverageExpression(data.sub.temp,group.by = "SID")$RNA %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% gather(cluster,ave_exp,-gene) %>% tbl_df() %>% mutate(species="mouse") %>% mutate(cluster=gsub("-","_",cluster))%>% rename(SID=cluster)
  
  
  saveRDS(species.miRNA.exp,file=paste0("tmp_data/",TD,"/","species.miRNA.ave.exp.rds"))
}

#' loading human anno
load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/Gene.meta.Rdata",verbose=T) 
human.miRNA.anno <- miRNA.bed %>% mutate(human=V7,miRNA=gsub("hsa-","",V7)) %>% select(human,miRNA) %>% unique()  %>% mutate(pmiRNA=sub("(miR-[0-9]+).*", "\\1", miRNA))
rm(mature_prec.anno,trans.anno,Gene.discrp,miRNA.family,miRNA.bed,chrX.miRNA,chrX.miRNA.unique,chrY.miRNA,chrY.miRNA.unique,chr.smallRNA.id)

load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)
mouse.miRNA.anno <- miRNA.bed %>% mutate(mouse=V7,miRNA=gsub("mmu-","",V7)) %>% select(mouse,miRNA) %>% unique()  %>% mutate(pmiRNA=sub("(miR-[0-9]+).*", "\\1", miRNA))
rm(mature_prec.anno,trans.anno,Gene.discrp,miRNA.family,miRNA.bed,chrX.miRNA,chrX.miRNA.unique,chr.smallRNA.id)


#' miRNA family 
mouse.miRNA.family <- readRDS(paste0("tmp_data/",TD,"/miRNA.family.rds")) %>% separate(miRNA_family,c("mf","ID"),sep=":") %>% select(-ID) %>% unique() 
human.miRNA.family <- read.delim("~/My_project/Extra_smncRNA/big_doc/hsa.mature_miRNA_family.out.tsv",stringsAsFactors = F,head=F) %>% tbl_df() %>% mutate(mf=paste(V2,V1,sep=":"),mature_miRNA=V3)%>% select(mf,mature_miRNA) %>% separate(mf,c("mf","ID"),sep=":") %>% select(-ID)  %>% unique()  

#' human and mouse share the same mature seq
hm.sameMature.out <- read.delim("~/Genome_new/mirBase/Release22.1/mouse.human.same.mature.fa.out",stringsAsFactors = F,sep="\t",head=F) %>% tbl_df() %>% select(V2,V4) %>% setNames(c("mouse","human")) %>% mutate(mouse_miRNA=gsub("mmu-","",mouse)) %>% unique()  %>% mutate(mouse_pmiRNA=sub("(miR-[0-9]+).*", "\\1", mouse_miRNA)) %>% mutate(human_miRNA=gsub("hsa-","",human)) %>% unique() %>% mutate(human_pmiRNA=sub("(miR-[0-9]+).*", "\\1", human_miRNA)) %>% unique() %>% filter(mouse %in% mouse.miRNA.anno$mouse) %>% filter(human %in% human.miRNA.anno$human)

hm.sameSeed.out <- read.delim("~/Genome_new/mirBase/Release22.1/mouse.human.same.seed.2_8p.fa.out",stringsAsFactors = F,sep="\t",head=F) %>% tbl_df() %>% select(V2,V3) %>% setNames(c("mouse","human")) %>% mutate(mouse_miRNA=gsub("mmu-","",mouse)) %>% unique()  %>% mutate(mouse_pmiRNA=sub("(miR-[0-9]+).*", "\\1", mouse_miRNA)) %>% mutate(human_miRNA=gsub("hsa-","",human)) %>% unique()  %>% mutate(human_pmiRNA=sub("(miR-[0-9]+).*", "\\1", human_miRNA))  %>% left_join(mouse.miRNA.family %>% mutate(mouse=mature_miRNA,mouse_mf=mf) %>% select(mouse,mouse_mf),by="mouse") %>% left_join(human.miRNA.family %>% mutate(human=mature_miRNA,human_mf=mf) %>% select(human,human_mf),by="human") %>% filter(mouse %in% mouse.miRNA.anno$mouse) %>% filter(human %in% human.miRNA.anno$human)

#' confident pairs ## same seed
conf.mf.pairs <- hm.sameSeed.out %>% filter(!is.na(mouse_mf) & !is.na(human_mf)) %>%select(mouse_mf,human_mf) %>% unique()
conf.pm.pairs <- hm.sameSeed.out  %>%select(mouse_pmiRNA,human_pmiRNA) %>% unique()
conf.mm.pairs <- hm.sameSeed.out  %>%select(mouse,human) %>% unique()


#' link mouse and human miRNA
#' by same mature sequence, 
hm.miRNA.anno <- hm.sameMature.out %>%  mutate(miRNA=ifelse(mouse_miRNA!=human_miRNA,paste(human_miRNA,mouse_miRNA,sep=":"),human_miRNA))  %>%  select(mouse,human,miRNA) %>% mutate(CF=3)  %>% unique()
nrow(hm.miRNA.anno )

#' by name and same seed sequence
hm.miRNA.anno<- mouse.miRNA.anno%>%select(-pmiRNA) %>% inner_join(human.miRNA.anno %>%select(-pmiRNA),by="miRNA",relationship = "many-to-many") %>% unique() %>% arrange(miRNA) %>% filter(!miRNA %in% hm.miRNA.anno$miRNA)  %>% inner_join(conf.mm.pairs,by =c("mouse", "human")) %>% mutate(CF=2) %>% bind_rows(hm.miRNA.anno)%>% unique()

#' choose the duplicated ones (filter multiple match)
hm.miRNA.anno <- hm.miRNA.anno %>% filter(! miRNA %in% c("miR-199b-3p:miR-199a-3p","miR-199a-3p:miR-199b-3p","miR-365a-3p:miR-365-3p","miR-365b-3p:miR-365-3p")) %>% filter(! (mouse=="mmu-miR-486b-5p" & human=="hsa-miR-486-5p"))

#' check duplicated
hm.miRNA.anno %>% group_by(mouse) %>% mutate(nC=n()) %>% filter(nC>1)
hm.miRNA.anno %>% group_by(human) %>% mutate(nC=n()) %>% filter(nC>1)


#' detail with family 
hm.miRNA.anno.detail <- hm.miRNA.anno%>% mutate(human_pmiRNA=sub("hsa-(miR-[0-9]+).*", "\\1", human),mouse_pmiRNA=sub("mmu-(miR-[0-9]+).*", "\\1", mouse)) %>% left_join(mouse.miRNA.family %>% mutate(mouse=mature_miRNA,mouse_mf=mf) %>% select(mouse,mouse_mf),by="mouse") %>% left_join(human.miRNA.family %>% mutate(human=mature_miRNA,human_mf=mf) %>% select(human,human_mf),by="human") %>% unique() 

#' fill the missing mf based on hm miRNA
hm.miRNA.anno.detail <- hm.miRNA.anno.detail  %>% mutate(mouse_mf=ifelse(is.na(mouse_mf),human_mf,mouse_mf)) %>% mutate(human_mf=ifelse(is.na(human_mf),mouse_mf,human_mf))
hm.miRNA.anno.detail %>% filter(mouse_mf!=human_mf)

mouse.miRNA.family <- hm.miRNA.anno.detail %>% filter(!is.na(mouse_mf)) %>% mutate(mf=mouse_mf,mature_miRNA=mouse) %>% select(mf,mature_miRNA) %>% unique() %>% bind_rows(mouse.miRNA.family) %>% unique()
human.miRNA.family <- hm.miRNA.anno.detail %>% filter(!is.na(human_mf)) %>% mutate(mf=human_mf,mature_miRNA=human) %>% select(mf,mature_miRNA) %>% unique() %>% bind_rows(human.miRNA.family) %>% unique()



#' link mf family
hm.mf.anno <- mouse.miRNA.family %>% mutate(mouse_mf=mf) %>% select(mouse_mf,mf) %>% unique() %>% inner_join(human.miRNA.family %>% mutate(human_mf=mf) %>% select(human_mf,mf) %>% unique(),by="mf") 
hm.mf.anno <- hm.mf.anno %>% inner_join(conf.mf.pairs,by = c("mouse_mf", "human_mf")) ## 26 family excluded
table(duplicated(hm.mf.anno $mouse_mf))
table(duplicated(hm.mf.anno $human_mf))

#' update family annotation
mouse.miRNA.family <- mouse.miRNA.family %>% mutate(mod_mf=ifelse(mf %in% hm.mf.anno$mf, paste0("ov_",mf), paste0("mouse_",mf)))

human.miRNA.family <- human.miRNA.family %>% mutate(mod_mf=ifelse(mf %in% hm.mf.anno$mf, paste0("ov_",mf), paste0("human_",mf)))


# split annotation
hm.miRNA.anno.stype <- mouse.miRNA.anno %>% rows_update(hm.miRNA.anno %>% select(miRNA,mouse) %>% filter(!is.na(mouse)),by="mouse" )%>% rename(ID=mouse) %>% left_join(mouse.miRNA.family %>% rename(ID=mature_miRNA),by="ID") %>% mutate(type=ifelse(ID %in% hm.miRNA.anno.detail$mouse,"ov","mouse_spec_miR")) %>% mutate(type=ifelse(type=="mouse_spec_miR" & mf %in% hm.mf.anno$mf ,"mouse_spec_member",type),species="mouse") %>% mutate(mod_miRNA=ifelse(type=="ov",paste0("ov_",miRNA),paste0("mouse_",miRNA))) %>% bind_rows(human.miRNA.anno %>% rows_update(hm.miRNA.anno %>% select(miRNA,human) %>% filter(!is.na(human)),by="human" )%>% rename(ID=human) %>% left_join(human.miRNA.family %>% rename(ID=mature_miRNA),by="ID") %>% mutate(type=ifelse(ID %in% hm.miRNA.anno.detail$human,"ov","human_spec_miR")) %>% mutate(type=ifelse(type=="human_spec_miR" & mf %in% hm.mf.anno$mf ,"human_spec_member",type),species="human")%>% mutate(mod_miRNA=ifelse(type=="ov",paste0("ov_",miRNA),paste0("human_",miRNA))) )

hm.miRNA.anno.stype %>% group_by(species,ID) %>% summarise(nT=n_distinct(type))  %>% filter(nT >1) #should be 0
table(hm.miRNA.anno.stype $type)


plot.results <- list()

#' check the embryo expressed miRNA family
temp.list <- list()
temp.list$mouse <- mouse.miRNA.family %>% filter(mature_miRNA %in% species.miRNA.exp$mouse.exp.miRNA) %>% pull(mod_mf) %>% unique()
temp.list$human <- human.miRNA.family %>% filter(mature_miRNA %in% species.miRNA.exp$human.exp.miRNA) %>% pull(mod_mf) %>% unique()


plot.results$hm.exp.miRNA.fov <- ggvenn::ggvenn(temp.list,fill_color=as.vector(species.col),show_percentage = F)+ggtitle("miRNA family")+FunTitle()
plot.results$hm.exp.miRNA.fov 


#' check the embryo expressed miRNA ID
temp.list <- list()

temp.list$mouse <- hm.miRNA.anno.stype %>% filter(ID %in% species.miRNA.exp$mouse.exp.miRNA) %>% filter(species=="mouse") %>% pull(mod_miRNA) %>% unique() 
temp.list$human <- hm.miRNA.anno.stype %>% filter(ID %in% species.miRNA.exp$human.exp.miRNA) %>% filter(species=="human") %>% pull(mod_miRNA) %>% unique() 
plot.results$hm.exp.pmiRNA.ov <- ggvenn::ggvenn(temp.list,fill_color=as.vector(species.col),show_percentage = F)+ggtitle("miRNA")+FunTitle()
plot.results$hm.exp.pmiRNA.ov

#' check the triangle marker expression in human and mouse
temp.sel.mf.anno <-  hm.miRNA.anno.stype  %>% filter(ID %in% species.miRNA.exp$human.t3.mk$gene | ID %in% species.miRNA.exp$mouse.t3.mk$gene) 

 
#' here requir those miRNA in each triangle are DE in one species
temp1 <- species.miRNA.exp$human.t3 %>% rename(ID=gene) %>% inner_join(temp.sel.mf.anno %>% filter(species=="human") %>% select(ID,miRNA,type),by="ID")%>% filter(SID %in% c("L8CM","ICM","TE")) %>% group_by(miRNA,SID,species,type) %>% summarise(Exp=sum(ave_exp)) 
temp2 <- species.miRNA.exp$mouse.t3 %>% rename(ID=gene) %>% inner_join(temp.sel.mf.anno%>% filter(species=="mouse") %>% select(ID,miRNA,type),by="ID") %>% filter(SID %in% c("L8CM","ICM","TE")) %>% group_by(miRNA,SID,species,type) %>% summarise(Exp=sum(ave_exp)) 
temp.exp <- temp1 %>% bind_rows(temp2) %>% ungroup() %>% rename(gene=miRNA) 

temp.hl.mf <- c("miR-127-3p","miR-370-3p","miR-409-3p","miR-369-3p","let-7f-5p","miR-27b-3p","miR-99b-5p","miR-378a-3p","miR-10b-5p","miR-140-3p","miR-652-3p","miR-146b-5p","miR-136-3p","miR-411-5p","miR-381-3p","miR-195-5p:miR-195a-5p","miR-200b-5p","miR-3116","miR-663b","miR-10396b-5p","miR-663a","miR-467a-3p","miR-466j","miR-669h-5p","miR-499a-5p:miR-499-5p","miR-466m-5p","miR-466b-3p","miR-10395-5p","miR-34a-3p","miR-181a-3p:miR-181a-1-3p")
#temp.pt.col <- c("human_spec_miR"="#e41a1c","human_spec_member"="#F8766D",ov="#00ba38","mouse_spec_miR"="#1f77b4","mouse_spec_member"="#00B0F6")
temp.pt.col <- c("human_spec_miR"="#e41a1c","human_spec_member"="#7A1FA2",ov="#00ba38","mouse_spec_miR"="#00B0F6","mouse_spec_member"="#CD9600")


temp.plot <- list()
for (n in c("human","mouse")) {
  temp.input <- temp.exp %>% filter(species==n)%>% group_by(gene,species) %>% mutate(sumExp=sum(Exp)) %>% mutate(propExp=Exp/sumExp) %>% select(gene,type,species,SID,propExp)  %>% spread(SID,propExp) %>% filter(!is.na(ICM))
  temp.hl.input <- temp.input %>% filter(gene %in% temp.hl.mf) %>% filter(ICM > 0.4 | TE >0.4 | L8CM > 0.4)
  temp.plot[[n]] <- ggtern(data=temp.input,aes(x=ICM,y=L8CM, z=TE)) + geom_point(mapping=aes(color=type),size=1.5)+annotate(geom='text',x=temp.hl.input$ICM,y=temp.hl.input$L8CM,z=temp.hl.input$TE, label =  temp.hl.input$gene ,color = temp.pt.col[temp.hl.input$type],size=2.5, vjust =c(1),alpha=0.75) +ggtern::theme_bw() + labs(title=n)+scale_color_manual(values=temp.pt.col)
  temp.plot[[n]]
}
plot.results$hm.miR.ggtern <- temp.plot
temp.plot$human
temp.plot$mouse



#' check dev related DEGs
load(paste0("tmp_data/",TD,"/lineage.segregation.DEG.out.Rdata"),verbose=T)
mouse.related.DE.miRNA  <-   DEG.results.list$miRNA[c("ICMs_vs_L8C","TEs_vs_L8C","L16C_vs_L8C","Early_ICM_vs_L16C","Late_ICM_vs_Early_ICM","EB_TE_vs_L16C" ,"MB_TE_vs_EB_TE","LB_TE_vs_MB_TE","LB_TE_vs_EB_TE","Early_ICM_vs_EB_TE","Early_ICM_vs_MB_TE","Late_ICM_vs_LB_TE","ICMs_vs_TEs")] %>% lapply(function(x) {x$DEG.result$gene}) %>% unlist() %>% unique()

load("~/My_project/Extra_smncRNA/tmp_data/Apr_2023/lineage.segregation.DEG.out.Rdata",verbose = T)
names(DEG.results)
human.related.DE.miRNA <- DEG.results  %>% lapply(function(x) {x$DEG.result$gene}) %>% unlist() %>% unique()


#' same ones
temp.sel.FM <- hm.miRNA.anno.stype %>% filter(type=="ov") %>% filter(ID %in% mouse.related.DE.miRNA | ID %in% human.related.DE.miRNA) %>% pull(miRNA) %>% unique() 

temp1 <- species.miRNA.exp$human.dev%>% rename(ID=gene) %>% inner_join(hm.miRNA.anno.stype ,by=c("species","ID")) %>% filter(miRNA %in% temp.sel.FM ) %>% group_by(miRNA,SID,species) %>% summarise(exp=sum(ave_exp)) %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"))
temp2 <- species.miRNA.exp$mouse.dev%>% rename(ID=gene) %>% inner_join(hm.miRNA.anno.stype ,by=c("species","ID")) %>% filter(miRNA %in% temp.sel.FM ) %>% group_by(miRNA,SID,species) %>% summarise(exp=sum(ave_exp)) %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"))

heat.col <- (colorRampPalette(c("royalblue3","white","firebrick4"))(50))
temp <- (temp1 %>%group_by(miRNA) %>% top_n(1,exp) %>% filter(exp >0) %>% inner_join(temp2 %>%group_by(miRNA) %>% top_n(1,exp) %>% filter(exp >0),by=c("miRNA","SID") ) %>% split(.,.$SID))[c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")] %>% do.call("bind_rows",.)
temp.plot <- list()
temp.plot$human.ph <- (temp1 %>% filter(miRNA %in% temp$miRNA) %>% select(-species) %>% spread(SID,exp) %>% select(miRNA,L8C,L16C,EB_ICM,MB_ICM,EB_TE,MB_TE) %>% tibble::column_to_rownames("miRNA") %>% log1p() %>% FunPreheatmapNoLog(Evalue=1.5))[temp$miRNA,] %>%t() %>% pheatmap(scale = "none",cluster_rows = F,cluster_cols = F,main="human",color=heat.col,cellwidth=10,cellheight=10 )%>%ggplotify::as.ggplot()
temp.plot$mouse.ph <- (temp2 %>% filter(miRNA %in% temp$miRNA) %>% select(-species) %>% spread(SID,exp) %>% select(miRNA,L8C,L16C,EB_ICM,MB_ICM,EB_TE,MB_TE) %>% tibble::column_to_rownames("miRNA") %>% log1p() %>% FunPreheatmapNoLog(Evalue=1.5))[temp$miRNA,] %>%t() %>% pheatmap(scale = "none",cluster_rows = F,cluster_cols = F,main="mouse",color=heat.col,cellwidth=10,cellheight=10 )%>%ggplotify::as.ggplot()

cowplot::plot_grid(temp.plot$human.ph,temp.plot$mouse.ph,ncol=1)
plot.results$hm.miR.same.ph <- temp.plot


#' diff pattern ones (species unique), 
temp1 <- species.miRNA.exp$human.dev %>% rename(ID=gene) %>% filter(ID %in% human.related.DE.miRNA) %>% inner_join(hm.miRNA.anno.stype ,by=c("species","ID")) %>% filter(type!="ov") %>% group_by(mf,type,miRNA,SID,species) %>% summarise(exp=sum(ave_exp)) %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")) 

temp2 <- species.miRNA.exp$mouse.dev%>% rename(ID=gene)%>% filter(ID %in% mouse.related.DE.miRNA) %>% inner_join(hm.miRNA.anno.stype ,by=c("species","ID")) %>% filter(type!="ov")  %>% group_by(mf,type,miRNA,SID,species) %>% summarise(exp=sum(ave_exp)) %>% filter(SID %in% c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE"))

temp.human <- (temp1 %>%group_by(miRNA,type,mf) %>% top_n(1,exp) %>% filter(exp >0) %>% anti_join(temp2 %>%group_by(miRNA,mf) %>% top_n(1,exp) %>% filter(exp >0) %>% select(mf,miRNA,SID),by=c("mf","miRNA","SID") ) %>% split(.,.$SID) %>% lapply(function(x){arrange(x,type,mf)}))[c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")] %>% do.call("bind_rows",.)%>% filter(!is.na(mf))
temp.mouse <- (temp2 %>%group_by(miRNA,type,mf) %>% top_n(1,exp) %>% filter(exp >0) %>% anti_join(temp1 %>%group_by(miRNA,mf) %>% top_n(1,exp) %>% filter(exp >0)%>% select(mf,miRNA,SID),by=c("mf","miRNA","SID") ) %>% split(.,.$SID)%>% lapply(function(x){arrange(x,type,mf)}))[c("L8C","L16C","EB_ICM","MB_ICM","EB_TE","MB_TE")] %>% do.call("bind_rows",.)%>% filter(!is.na(mf))



temp.sel.FM.family <- temp.human %>% bind_rows(temp.mouse) %>% filter(!is.na(mf))  %>% group_by(mf,SID,species) %>% summarise(nM=n_distinct(miRNA),exp=sum(exp)) %>% filter(nM >0) # all

temp.sel.FM.family.label <- temp.sel.FM.family  %>% mutate(SSID=paste(species,SID,sep=":")) %>% split(.,.$SSID) %>% lapply(function(x) {x %>% arrange(desc(nM),desc(exp)) %>% head(3)}) %>% do.call("bind_rows",.) %>% mutate(mf=paste0(mf," (",nM,")")) %>% filter(nM >1)%>% group_by(SSID) %>% summarise(mf=paste(mf,collapse=","))
temp.sel.FM.family.label

temp.sel.FM.family.hl <- temp.sel.FM.family  %>% mutate(SSID=paste(species,SID,sep=":")) %>% split(.,.$SSID) %>% lapply(function(x) {x %>% arrange(desc(nM),desc(exp)) %>% head(3)}) %>% do.call("bind_rows",.) %>% select(mf,species) %>% unique()
temp.sel.FM.family.hl.human <- temp.sel.FM.family.hl   %>% filter(species=="human") %>% pull(mf) %>% unique()
temp.sel.FM.family.hl.mouse <- temp.sel.FM.family.hl   %>% filter(species=="mouse") %>% pull(mf) %>% unique()


temp.human <- temp.human %>% inner_join(temp.sel.FM.family %>%select(mf,SID,species),by=c("mf","SID","species")) # no use
temp.mouse <- temp.mouse %>% inner_join(temp.sel.FM.family %>%select(mf,SID,species),by=c("mf","SID","species")) # no use

temp.human.mm.anno <- temp.human %>% ungroup()%>% select(miRNA,type,mf) %>% unique() %>% mutate(mf=ifelse(mf %in% temp.sel.FM.family.hl.human ,mf,"other")) %>% tibble::column_to_rownames("miRNA")
temp.mouse.mm.anno <- temp.mouse%>% ungroup() %>% select(miRNA,type,mf) %>% unique() %>% mutate(mf=ifelse(mf %in% temp.sel.FM.family.hl.mouse ,mf,"other")) %>% tibble::column_to_rownames("miRNA")


##### need to confirm the mf has the color
unique(c(temp.sel.FM.family.hl.human ,temp.sel.FM.family.hl.mouse))
temp.sel.FM.type.color.anno <- list(type=temp.pt.col[c("human_spec_miR","human_spec_member","mouse_spec_miR","mouse_spec_member")],mf=c("mir-515"="#F8766D","mir-290"="#EC8239","mir-506"="#DB8E00","mir-373"="#C79800","mir-1268"="#AEA200","mir-3180"="#8FAA00","mir-500"="#64B200","mir-743"="#00B81B","mir-548"="#00BD5C","mir-302"="#00C085","mir-28"="#00C1A7","mir-1307"="#00BFC4","mir-467"="#00BADE","mir-1193"="#00B2F3","mir-466"="#00A6FF","mir-298"="#7C96FF","mir-1843"="#B385FF","mir-7"="#D874FD","mir-465"="#EF67EB","mir-883"="#FD61D3","mir-672"="#FF63B6","mir-351"="#FF6B94","other"="grey77"))

temp.plot <- list()
temp.plot$human.ph <- (temp1 %>% filter(miRNA %in% temp.human$miRNA)  %>% ungroup()%>% select(-c(species,mf,type)) %>% spread(SID,exp) %>% select(miRNA,L8C,L16C,EB_ICM,MB_ICM,EB_TE,MB_TE) %>% tibble::column_to_rownames("miRNA") %>% log1p() %>% FunPreheatmapNoLog(Evalue=1.5))[temp.human$miRNA,] %>%t() %>% pheatmap(scale = "none",cluster_rows = F,cluster_cols = F,main="human(specie-spec)",color=heat.col,show_colnames=F,annotation_col=temp.human.mm.anno,annotation_colors=temp.sel.FM.type.color.anno)%>%ggplotify::as.ggplot()# ,,cellwidth=10,cellheight=10
temp.plot$mouse.ph <- (temp2 %>% filter(miRNA %in% temp.mouse$miRNA)  %>% ungroup()%>% select(-c(species,mf,type)) %>% spread(SID,exp) %>% select(miRNA,L8C,L16C,EB_ICM,MB_ICM,EB_TE,MB_TE) %>% tibble::column_to_rownames("miRNA") %>% log1p() %>% FunPreheatmapNoLog(Evalue=1.5))[temp.mouse$miRNA,] %>%t() %>% pheatmap(scale = "none",cluster_rows = F,cluster_cols = F,main="mouse(specie-spec)",color=heat.col,show_colnames=F ,annotation_col=temp.mouse.mm.anno,annotation_colors=temp.sel.FM.type.color.anno)%>%ggplotify::as.ggplot()#annotation_col=temp.mouse.mf.anno,
#,cellwidth=10,cellheight=10
cowplot::plot_grid(temp.plot$human.ph,temp.plot$mouse.ph,ncol=1)
plot.results$hm.miR.diff.species.ph <- temp.plot
