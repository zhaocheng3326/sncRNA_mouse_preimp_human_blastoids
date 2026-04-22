#' ---
#' title: traj inferring for embryonic cells 
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
  library(monocle) 
  #library(batchelor)
  #library(SeuratWrappers)
  library(slingshot)
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

load(paste0("tmp_data/",TD,"/Gene.meta.Rdata"),verbose = T)


sel.type <- c("miRNA")
meta.filter <- readRDS(paste0("tmp_data/",TD,"/small.meta.filter.rds"))  %>% mutate(pj=batch) 
meta.filter <- meta.filter %>% filter(batch !="batch4") %>% bind_rows(meta.filter %>% filter(batch=="batch4" & stage %in% batch4.sel.stage))

counts.filter <- readRDS(paste0("tmp_data/",TD,"/small.counts.filter.rds"))
miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(counts.filter$mature$ID)
chr.miRNA.ID <- trans.anno$mature %>% filter(type=="miRNA") %>% pull(ID) %>% intersect(chr.smallRNA.id) %>% intersect(counts.filter$mature$ID)
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% chr.miRNA.ID) %>% tibble::column_to_rownames("ID"))[,meta.filter$cell]

#' update the coseq small part annotation
meta.filter <- meta.filter %>% mutate(EML=devTime,RNA_EML="None",small_EML="None") %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Msmall.coseq_small.data.ob.umap.rds")) %>% select(cell,EML,RNA_EML,small_EML),by="cell")

#' update the full cell annotation 
meta.filter <- meta.filter  %>% rows_update(readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".data.ob.umap.rds")) %>% select(cell,EML),by="cell")

#' running
sel.type <- c("miRNA")

#'miRNA only
counts.miRNA.filter <- (counts.filter$mature %>% filter(ID %in% miRNA.ID) %>% tibble::column_to_rownames("ID")) [,meta.filter$cell]


#' normalization values using all-cell norm values
sel.exp <- readRDS(paste0("tmp_data/",TD,"/","allCells.miRNA",".norm.rds"))

data.ob <- readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.rds"))
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/","main.miRNA",".data.ob.umap.rds")) %>% inner_join(meta.filter,by="cell")

#' traj analysis
if (file.exists(paste0("tmp_data/",TD,"/mouse.miRNA.SCP.main.traj.Rdata"))) {
  load(paste0("tmp_data/",TD,"/mouse.miRNA.SCP.main.traj.Rdata"),verbose=T)
}else{
  data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(EML!="unknown") %>% pull(cell)))
  data.temp@meta.data$SID <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"]
  
  data.sce <- as.SingleCellExperiment(data.temp)
  
  data.sce <- slingshot(data.sce,  reducedDim = 'UMAP',clusterLabels = 'EML', start.clus = 'L2and4C',end.clus=c("ICM","TE"))
  data.crv <- SlingshotDataSet(data.sce)
  plot(reducedDims(data.sce)$UMAP)
  lines(data.crv, lwd=2, type = 'lineages', col = 'black')
  plot(reducedDims(data.sce)$UMAP,col="grey",pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab = "", ylab = "", axes = FALSE)
  lines(data.crv, lwd=2, type = 'lineages', col = c("royalblue"))
  
  
  data.pseudotime.na <-  slingPseudotime(data.sce, na=T) %>% as.data.frame() %>% setNames(c("ICM_traj","TE_traj")) %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% inner_join(data.ob.umap %>% select(cell,devTime,stage,EML),by="cell")
  data.pseudotime <-  slingPseudotime(data.crv, na=FALSE) %>% as.data.frame()  %>% setNames(c("ICM_traj","TE_traj"))%>% tibble::rownames_to_column("cell") %>% tbl_df() %>%  gather(traj,psdt,-cell) %>% inner_join(data.ob.umap %>% select(cell,devTime,stage,EML),by="cell")
  data.pseudotime.ave <- slingAvgPseudotime(data.crv)
  data.cellWeights <- slingCurveWeights(data.crv)
  
  data.pseudotime.mod <- as.data.frame(data.pseudotime.ave) %>% tibble::rownames_to_column("cell") %>% setNames(c("cell","ICM_traj")) %>% tbl_df() %>% mutate(TE_traj=ICM_traj) %>%  gather(traj,psdt,-cell) %>% inner_join(data.ob.umap %>% select(cell,devTime,stage,EML),by="cell")
  
  #data.pseudotime.ave.mod <- as.data.frame(data.pseudotime.ave) %>%  tibble::rownames_to_column("cell")  %>% tbl_df() %>% setNames(c("cell","PSDT")) %>% inner_join(data.ob.umap %>% select(cell,EML),by="cell")
  #data.pseudotime.mod <- data.pseudotime.mod 
  
  #'scale the pseudo time to align the time
  temp <- as.data.frame(data.pseudotime.ave) %>%  tibble::rownames_to_column("cell")  %>% tbl_df() %>% setNames(c("cell","PSDT")) %>% inner_join(data.ob.umap %>% select(cell,EML),by="cell")
  temp.cut <- temp %>% filter(EML=="L8CM") %>% arrange(PSDT) %>% tail(1) %>% pull(PSDT)
  temp.TE <- temp %>% filter(EML=="TE" & PSDT > temp.cut) %>% mutate(PSDT=PSDT-temp.cut) %>% inner_join( data.pseudotime.mod %>% select(cell,devTime) %>% unique(),by="cell")
  temp.ICM <- temp %>% filter(EML=="ICM" & PSDT > temp.cut) %>% mutate(PSDT=PSDT-temp.cut)%>% inner_join( data.pseudotime.mod %>% select(cell,devTime)%>% unique(),by="cell")
  #temp.scale.ft <- max(temp.ICM$PSDT)/max(temp.TE$PSDT)
  temp.scale.ft <- mean(temp.ICM %>% filter(devTime=="64C") %>% pull(PSDT))/mean(temp.TE %>% filter(devTime=="64C") %>% pull(PSDT))
  print(temp.scale.ft)
  temp.TE <- temp.TE %>% mutate(PSDT=PSDT * temp.scale.ft+temp.cut)
  data.pseudotime.ave.mod <- temp %>% rows_update(temp.TE %>% select(-devTime),by="cell")
  data.pseudotime.mod <- data.pseudotime.mod %>% rows_update(temp.TE %>% mutate(psdt=PSDT) %>% select(cell,psdt),by="cell")
  
  
  #' select EML
  traj_lineage_sel <- data.frame(EML=c("L2and4C","L8CM","ICM"),traj="ICM_traj",sel_traj="ICM_main_traj")  %>% bind_rows(data.frame(EML=c("L2and4C","L8CM","TE"),traj="TE_traj",sel_traj="TE_main_traj")) %>% tbl_df() %>% mutate_all(as.vector)
  
  data.pseudotime.sel <- data.pseudotime.mod %>% inner_join(traj_lineage_sel,by=c("traj","EML"),multiple = "all")
  
  #traj_lineage_sel <- data.frame(EML=c("oocyte","L2and4C","L8CM","ICM"),traj="ICM_traj",sel_traj="ICM_main_traj")  %>% bind_rows(data.frame(EML=c("oocyte","L2and4C","L8CM","TE"),traj="TE_traj",sel_traj="TE_main_traj")) %>% tbl_df() %>% mutate_all(as.vector)
  #data.pseudotime.sel <- data.pseudotime.mod %>% inner_join(traj_lineage_sel,by=c("traj","EML"),multiple = "all")
  #' modified from SCP tools https://github.com/zhanghao-njmu/SCP
  data.temp <- subset(data.ob,cell=unique(data.pseudotime.sel$cell))
  data.temp@meta.data <- data.temp@meta.data %>% cbind((data.pseudotime.sel %>% select(cell,psdt,sel_traj) %>% spread(sel_traj,psdt) %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),])
  psdt.genes <- list()
  res.list <- list()
  # modified from SCP pacakage
  for (stj in unique(data.pseudotime.sel$sel_traj)) {
    temp.pseudotime <- data.pseudotime.sel %>% filter(sel_traj==stj)  %>% select(cell,psdt,sel_traj,EML,devTime) %>% arrange(psdt)
    temp.pseudotime.vec <- temp.pseudotime %>% select(cell,psdt) %>% tibble::column_to_rownames("cell") %>% as.matrix()
    temp.sel.exp <- sel.exp[,temp.pseudotime$cell]
    temp.sel.exp <- temp.sel.exp[rowSums(temp.sel.exp >0 ) >= 5,]
    gam_out <- list()
    Y_ordered <- temp.sel.exp %>% as.matrix()
    t_ordered <- temp.pseudotime.vec[,"psdt"]
    for (n in seq_len(nrow(Y_ordered))) {
      print(paste0(stj,n))
      feature_nm <- rownames(Y_ordered)[n]
      family_use <- "gaussian"
      sizefactror <- 1
      
      mod <- mgcv::gam(y ~ s(x, bs = "cs") + offset(rep(log(1),ncol(Y_ordered))),family = rep(family_use,ncol(Y_ordered)),data = data.frame(y = Y_ordered[feature_nm,,drop=T ], x = t_ordered))
      
      pre <- predict(mod, type = "link", se.fit = TRUE)
      upr <- pre$fit + (2 * pre$se.fit)
      lwr <- pre$fit - (2 * pre$se.fit)
      upr <- mod$family$linkinv(upr)
      lwr <- mod$family$linkinv(lwr)
      res <- summary(mod)
      fitted <- fitted(mod)
      pvalue <- res$s.table[[4]]
      dev.expl <- res$dev.expl
      r.sq <- res$r.sq
      fitted.values <- fitted * sizefactror
      upr.values <- upr * sizefactror
      lwr.values <- lwr * sizefactror
      exp_ncells <- sum(Y_ordered[feature_nm, ] > min(Y_ordered[feature_nm, ]), na.rm = TRUE)
      peaktime <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
      valleytime <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])
      gam_out[[feature_nm]] <- list(
        features = feature_nm, exp_ncells = exp_ncells,
        r.sq = r.sq, dev.expl = dev.expl,
        peaktime = peaktime, valleytime = valleytime,
        pvalue = pvalue, fitted.values = fitted.values,
        upr.values = upr.values, lwr.values = lwr.values
      )
    }
    
    raw_matrix <- Y_ordered
    fitted_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["fitted.values"]]))
    colnames(fitted_matrix) <- rownames(Y_ordered)
    fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)
    
    upr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["upr.values"]]))
    colnames(upr_matrix) <- rownames(Y_ordered)
    upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)
    
    lwr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["lwr.values"]]))
    colnames(lwr_matrix) <- rownames(Y_ordered)
    lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)
    DynamicFeatures <- as.data.frame(do.call(rbind.data.frame, lapply(gam_out, function(x) x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")]))) %>% filter(!is.na(peaktime))
    char_var <- c("features")
    numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
    DynamicFeatures[, char_var] <- lapply(DynamicFeatures[, char_var, drop = FALSE], as.character)
    DynamicFeatures[, numb_var] <- lapply(DynamicFeatures[, numb_var, drop = FALSE], as.numeric)
    rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
    DynamicFeatures[, "padjust"] <- p.adjust(DynamicFeatures[, "pvalue", drop = TRUE])
    raw_matrix <- raw_matrix[rownames(DynamicFeatures),]
    fitted_matrix <- fitted_matrix[,c("pseudotime",rownames(DynamicFeatures))]
    upr_matrix <- upr_matrix[,colnames(fitted_matrix)]
    lwr_matrix <- upr_matrix[,colnames(fitted_matrix)]
    
    res <- list(
      DynamicFeatures = DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = fitted_matrix,
      upr_matrix = upr_matrix,
      lwr_matrix = lwr_matrix,
      libsize = NULL,
      lineages = stj,
      family = family_use
    )
    data.temp@tools[[paste0("DynamicFeatures_", stj)]] <- res
    
    temp.res <- res$DynamicFeatures %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% filter(padjust < 0.05 & exp_ncells > 10 & r.sq > 0.2 & dev.expl > 0.2) %>% arrange(peaktime) ## cutoff from scp
    k <- 6;set.seed(123);temp.ph <- temp.sel.exp[temp.res$gene,temp.pseudotime$cell] %>% FunPreheatmapNoLog() %>% pheatmap::pheatmap(scale="none",cluster_rows = F,cluster_cols = F,kmeans_k = k,show_rownames = F,show_colnames = F)
    temp.res <- temp.res %>% inner_join(as.data.frame(temp.ph$kmeans$cluster) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% setNames(c("gene","kcluster")) %>% mutate(kcluster=paste0("K",kcluster)),by="gene")
    temp.sel.exp[temp.res$gene,temp.pseudotime$cell] %>% FunPreheatmapNoLog() %>% pheatmap::pheatmap(scale="none",cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_rows=(temp.res %>% select(gene,cluster) %>% tibble::column_to_rownames("gene")),annotation_col = (temp.pseudotime %>% tibble::column_to_rownames("cell")),main=stj,annotation_row = (temp.res %>% select(gene,cluster) %>% tibble::column_to_rownames("gene")))
  
    psdt.genes[[stj]] <- temp.res
    res.list[[stj]] <- res
  }
  
  save(psdt.genes,data.sce,data.crv,data.pseudotime.sel, data.pseudotime.na,data.pseudotime.mod,res.list,file=paste0("tmp_data/",TD,"/mouse.miRNA.SCP.main.traj.Rdata") )
}
