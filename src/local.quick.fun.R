EML.od <- c("sperm","oocyte","L2and4C","L8CM","TE","ICM")
devTime.od <- c("sperm","oocyte","2C","4C","8C","16C","32C","64C")
batch4.sel.stage <- c("oocyte","0.1 uL sperm aliqout","16 cell","32 cell blast")

#batch5.lq.embyro.ID <- c("em_LTT13","em_LTT14","em_LTT7","em_LTV10","em_LTV12","em_LTV4","em_YTT74","em_YTV72","em_YTV73")
split2.lq.embyro.ID <- c("em_Q9","em_R1","em_R2","em_R3","em_R4")
lq.small.cells <- c("Sample_1727_0")

APPJ_UMAP <- function(temp.umap,PJ) {
  temp.plot <- list()
  for ( n in unique(temp.umap %>% filter(pj==PJ) %>% pull(EML))) {
    temp.plot[[n]] <-ggplot()+ geom_point(temp.umap %>% filter(!(pj==PJ & EML==n)),mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.25)+geom_point( temp.umap %>% filter((pj==PJ & EML==n)),,mapping=aes(x=UMAP_1,y=UMAP_2),col="red",size=0.25 )+ggtitle(n)+theme_void()
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}
APPJ_devTime_UMAP <- function(temp.umap,PJ) {
  temp.plot <- list()
  for ( n in unique(temp.umap %>% filter(pj==PJ) %>% pull(devTime))) {
    temp.plot[[n]] <-ggplot()+ geom_point(temp.umap %>% filter(!(pj==PJ & devTime==n)),mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.25)+geom_point( temp.umap %>% filter((pj==PJ & devTime==n)),,mapping=aes(x=UMAP_1,y=UMAP_2),col="red",size=0.25 )+ggtitle(n)+theme_void()
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}
AP_UMAP <-  function(data.temp.UMAP) {
  temp.plot <-list()
  temp.plot$EML <-   ggplot()+geom_point(data.temp.UMAP ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.5)+geom_text(data=data.temp.UMAP %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=3 ) + theme_classic()+NoAxes()+NoAxes()+ggtitle("EML")+theme(plot.title = element_text(hjust=0.5))+NoLegend()
  #temp.plot$devTime <-   ggplot()+geom_point(data.temp.UMAP ,mapping=aes(x=UMAP_1,y=UMAP_2,color=devTime),size=0.5)+geom_text(data=data.temp.UMAP %>% group_by(devTime) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=devTime),size=3 ) + theme_classic()+NoAxes()+NoAxes()+ggtitle("devTime")+theme(plot.title = element_text(hjust=0.5))+NoLegend()
  print(cowplot::plot_grid(plotlist=temp.plot))
}
APPJ <- function(data.temp,PJ) {
  temp.plot <- list()
  for ( n in unique(data.temp@meta.data$EML[data.temp@meta.data$pj==PJ])) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj==PJ & data.temp@meta.data$EML==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}
APPJ_devTime <- function(data.temp,PJ) {
  temp.plot <- list()
  for ( n in unique(data.temp@meta.data$devTime[data.temp@meta.data$pj==PJ])) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj==PJ & data.temp@meta.data$devTime==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}

FunDEG_custom <- function(data.temp,cells.list,temp.compair) {
  s1=unlist(strsplit(temp.compair,split="_vs_"))[1]
  s2=unlist(strsplit(temp.compair,split="_vs_"))[2]
  c1=cells.list[[s1]]
  c2=cells.list[[s2]]
  data.temp <- subset(data.temp,cells=c(c1,c2))
  data.temp@meta.data[c(c1,c2),"id"]=c(rep("c1",length(c1)), rep("c2",length(c2)))
  Idents(data.temp) <- factor(data.temp@meta.data$id)
  DefaultAssay(data.temp) <- "RNA"
  
  G1G2.DEG <- FindMarkers(data.temp,ident.1="c1",ident.2="c2",verbose = F,test.use="wilcox",logfc.threshold=0.1) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% mutate(fdr=p.adjust(p_val,method="BH",n=nrow(data.temp)))
  #G1G2.DEG.sig <- G1G2.DEG %>% filter(p_val_adj < 0.05)
  #G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0 & p_val_adj <0.05)
  #G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0 & p_val_adj <0.05)
  
  G1G2.DEG.sig <- G1G2.DEG %>% filter(fdr < 0.05 & abs(avg_log2FC) >0.1)
  G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0.1 & fdr <0.05)
  G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0.1 & fdr <0.05)
  
  #G1G2.DEG.sig <- G1G2.DEG %>% filter(fdr < 0.05 & abs(avg_log2FC) >0.25)
  #G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0.25 & fdr <0.05)
  #G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0.25 & fdr <0.05)
  
  temp.out <- list()
  temp.out[["DEG.all.result"]] <- G1G2.DEG
  temp.out[["DEG.result"]] <- G1G2.DEG.sig
  temp.out[["DEG.result.up"]] <-  G1.sig.up
  temp.out[["DEG.result.down"]] <-  G2.sig.up
  return(temp.out)
}
APHL_UMAP <-  function(data.temp.UMAP,temp.sel.cells) {
  ggplot()+geom_point(data.temp.UMAP %>% filter(!cell %in% temp.sel.cells),mapping=aes(x=UMAP_1,y=UMAP_2),color="grey")+geom_point(data.temp.UMAP %>% filter(cell %in% temp.sel.cells),mapping=aes(x=UMAP_1,y=UMAP_2),color="red")
}
FunDEG_custom_strict <- function(data.temp,cells.list,temp.compair) {
  s1=unlist(strsplit(temp.compair,split="_vs_"))[1]
  s2=unlist(strsplit(temp.compair,split="_vs_"))[2]
  c1=cells.list[[s1]]
  c2=cells.list[[s2]]
  data.temp <- subset(data.temp,cells=c(c1,c2))
  data.temp@meta.data[c(c1,c2),"id"]=c(rep("c1",length(c1)), rep("c2",length(c2)))
  Idents(data.temp) <- factor(data.temp@meta.data$id)
  DefaultAssay(data.temp) <- "RNA"
  
  G1G2.DEG <- FindMarkers(data.temp,ident.1="c1",ident.2="c2",verbose = F,test.use="wilcox",logfc.threshold=0.1) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% mutate(fdr=p.adjust(p_val,method="BH",n=nrow(data.temp)))
  
  G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0.25 & fdr <0.05 & pct.1 > 0.5 & pct.2 < 0.1)
  G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0.25 & fdr <0.05 & pct.2 > 0.5 & pct.1 < 0.1 )
  G1G2.DEG.sig <- G1.sig.up %>% bind_rows(G2.sig.up)
  
  temp.out <- list()
  temp.out[["DEG.all.result"]] <- G1G2.DEG
  temp.out[["DEG.result"]] <- G1G2.DEG.sig
  temp.out[["DEG.result.up"]] <-  G1.sig.up
  temp.out[["DEG.result.down"]] <-  G2.sig.up
  return(temp.out)
}

FunTitle <- function() {
  p <- theme(plot.title = element_text(hjust=0.5))
  return(p)
}

#co.mk <- c("mmu-miR-99a-5p","mmu-miR-320-3p","mmu-miR-27a-3p","mmu-let-7b-5p","mmu-miR-22-3p","mmu-miR-23b-3p","mmu-miR-23a-3p","mmu-miR-221-3p", "mmu-miR-205-3p","mmu-miR-203-3p","mmu-miR-144-3p","mmu-miR-181a-5p","mmu-miR-205-5p","mmu-miR-130a-3p")
co.mk <- c("mmu-miR-99a-5p","mmu-miR-320-3p","mmu-miR-27a-3p","mmu-let-7b-5p","mmu-miR-22-3p","mmu-miR-23b-3p","mmu-miR-23a-3p","mmu-miR-221-3p", "mmu-miR-203-3p","mmu-miR-181a-5p","mmu-miR-205-5p","mmu-miR-130a-3p") %>% unique() #remove "mmu-miR-205-3p","mmu-miR-21a-5p","mmu-miR-27b-3p","mmu-miR-144-3p",

tp.mk <- c("mmu-miR-541-5p","mmu-miR-409-3p","mmu-miR-378a-3p","mmu-miR-10b-5p","mmu-miR-871-3p","mmu-miR-140-3p")

w.mk <- c("mmu-miR-302d-5p","mmu-miR-302d-3p","mmu-miR-302b-3p")


#small.lineage.col.set <- c("Loocyte"="#FCC1D7","L2C"="#F8766D", L4C="#00BF7D",L8C="#00CF7D",L16C="#A3A500",early32C="#7CAE00",early_ICM="#FB61D7",late_ICM="#A58AFF", early_TE="#00C6EB", late_TE ="#00B0F6",unknown="grey66")




#lineage.col.set <- c("sperm"="#00BE67","oocyte"="#F8766D",L2and4C="#CD9600",L8CM="#00BFC4", Prelineage="#C49A00",ICM="#FB61D7",TE="#00B0F6","L8CM"="#C49A00","L2and4C"="#00BF7D","EPI"="#FFAACF",PE="#7A1FA2","EB_TE"="#80D8FF","MB_TE"="#00CBF6","LB_TE"="#005A7B","unknown"="grey66")

lineage.col.set <- c("sperm"="#00BE67","oocyte"="#F8766D",L2and4C="#CD9600",L8CM="#00BFC4", Prelineage="#C49A00",ICM="#FB61D7",TE="#00B0F6","L8CM"="#C49A00","L2and4C"="#00BF7D","Early_ICM"="#FFAACF","Late_ICM"="#7A1FA2","EB_TE"="#B3E5E1","MB_TE"="#00CBF9","LB_TE"="#005A7B","unknown"="grey66","HTE"="#0085EE")


type.col.set <- c(miRNA="#F8766D",piRNA="#B79F00",rRNA="#00BA38",snoRNA="#00BFC4",snRNA="#619CFF",tRNA="#F564E3")

devTime.col.set <- c("sperm"="#00BE67","oocyte"="#F8766D","2C"="#CD9600","4C"="#7CAE00","8C"="#00BFC4","16C"="#C77CFF","32C"="#00A9FF","64C"="#FF61CC")
devTime.col.set2 <- c("Zygote"="#F8766D","2C"="#CD9600","4C"="#7CAE00","8C"="#00BFC4","16C"="#C77CFF","EB"="#00BE67","MB"="#00A9FF","LB"="#FF61CC")


devTime.col.set3 <- c("oocyte_4C"=as.vector(devTime.col.set["oocyte"]),"8C_morula"="#00BFC4",lineage.col.set [c("ICM","TE")])

devTime.od <- c("sperm","oocyte","2C","4C","8C","16C","32C","64C")
devTime.od2 <- c("oocyte","2C","4C","8C","16C","32C","64C")

#heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
heat.col <- colorRampPalette(c("#09005EFF","#0D0887FF","#0D0887FF","#110788","#1C078B","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
heat.col2 <- colorRampPalette(c("#09005EFF","#0D0887FF","#0D0887FF","#110788","#1C078B","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF","#F0F921FF"))(100)

species.col <- c("mouse"="#1F77B4","human"="#e41a1c")
Group.col <- c("TT.Lats1_2"="#00B0E1","TT.Yap_Taz"="#FB61A1", "TV.Lats1_2"="#4DAF4A","TV.Yap_Taz"="#A6D854")


HE.HB.NP.col.set <- c(lineage.col.set[c("ICM","TE","Early_ICM","MB_TE")] %>% setNames(c("ICM","TE","ELC","TLC")),"Naive"="#00CC55","Primed"="#E26EF7")


Fun_ModuleRadarPlot_mod <- function (seurat_obj, group.by = NULL, barcodes = NULL, combine = TRUE, 
                                     ncol = 4, wgcna_name = NULL, fill = TRUE, draw.points = FALSE, 
                                     ...) 
{
  
  seurat_obj <- data.query.perum
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  if (!require("ggradar")) {
    print("Missing package: ggradar")
    print("Installing package: ggradar")
    devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
  }
  meta <- seurat_obj@meta.data
  if (is.null(group.by)) {
    cell_grouping <- Idents(seurat_obj)
  }
  else {
    cell_grouping <- seurat_obj@meta.data[, group.by]
    names(cell_grouping) <- colnames(seurat_obj)
  }
  if (is.factor(cell_grouping)) {
    group_order <- levels(cell_grouping)
  }
  else {
    group_order <- unique(cell_grouping)
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  mod_colors <- modules %>% dplyr::select(c(module, color)) %>% 
    dplyr::distinct()
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  MEs <- GetMEs(seurat_obj)
  MEs <- MEs[, colnames(MEs) != "grey"]
  MEs <- MEs %>% as.data.frame()
  if (!is.null(barcodes)) {
    if (!(all(barcodes %in% colnames(seurat_obj)))) {
      stop("Invalid selection for barcodes, some are not found in the colnames(seurat_obj)")
    }
    MEs <- MEs[barcodes, ]
    meta <- meta[barcodes, ]
    cell_grouping <- cell_grouping[barcodes]
  }
  MEs$cluster <-cell_grouping
  clusters <- as.character(unique(cell_grouping))
  plot_df <- MEs %>% dplyr::group_by(cluster) %>% dplyr::summarise_all(mean) %>% 
    as.data.frame()
  rownames(plot_df) <- plot_df$cluster
  plot_df <- dplyr::select(plot_df, -cluster)
  plot_df <- t(plot_df) %>% as.data.frame()
  plot_df[plot_df < 0] <- 0
  plot_df$group <- rownames(plot_df)
  plot_df <- plot_df[, c("group", clusters)]
  print(head(plot_df))
  plot_df$group <- factor(as.character(plot_df$group), levels = mods)
  plot_df <- plot_df %>% dplyr::arrange(group) %>% as.data.frame()
  colnames(plot_df) <- c("group", clusters)
  plot_df <- plot_df[, c("group", group_order)]
  plot_list <- list()
  for (i in 1:nrow(plot_df)) {
    cur_mod <- as.character(plot_df[i, "group"])
    cur_color <- subset(mod_colors, module == cur_mod) %>% 
      .$color
    plot_list[[cur_mod]] <- ggradar::ggradar(plot_df[i, ], 
                                             group.colours = cur_color, draw.points = draw.points, 
                                             fill = fill) + Seurat::NoLegend() + ggtitle(cur_mod) + 
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  if (combine) {
    patch <- patchwork::wrap_plots(plot_list, ncol)
    return(patch)
  }
  else {
    return(plot_list)
  }
}


FunMiRNA_TF_fisher <- function(sel_miRNA,TF_miRNA) {
  bg_miRNA <- TF_miRNA$miRNA %>% unique() 
  sel_miRNA <- sel_miRNA %>% intersect(bg_miRNA )%>% unique()
  
  temp <- TF_miRNA %>% filter(miRNA %in% sel_miRNA) %>% group_by(TF) %>% summarise(TSet=n_distinct(miRNA))  %>% mutate(BgSet=length(sel_miRNA)) %>% inner_join(TF_miRNA  %>% group_by(TF) %>% summarise(TBg=n_distinct(miRNA),Bg=length(bg_miRNA)) ,by="TF") 
  
  temp$stat <- temp %>% apply(1,function(x){
    temp.fish.cal <- fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE),alternative = "great")
    if (temp.fish.cal$estimate=="Inf") {
      temp.fish.cal.mod <- fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE)+1,alternative = "great")
      return(paste(temp.fish.cal.mod$estimate,temp.fish.cal$p.value,sep=":"))
    }else{
      return(paste(temp.fish.cal$estimate,temp.fish.cal$p.value,sep=":"))
    }
    
  })
  temp <- temp %>% separate(stat,c("od","pvalue"),sep = ":") %>% mutate(od=as.numeric(od),pvalue=as.numeric(pvalue)) %>% mutate(p_val_adj=p.adjust(pvalue,method="fdr"))
  return(temp)

}

FunGeneSet_miRNA_fisher <- function(sel_gene,miRNA_TG) {
  bg_gene <- miRNA_TG$gene %>% unique() 
  sel_gene <- sel_gene %>% intersect(bg_gene )%>% unique()
  
  temp <- miRNA_TG %>% filter(gene %in% sel_gene) %>% group_by(miRNA) %>% summarise(TSet=n_distinct(gene))  %>% mutate(BgSet=length(sel_gene)) %>% inner_join(miRNA_TG  %>% group_by(miRNA) %>% summarise(TBg=n_distinct(gene),Bg=length(bg_gene)) ,by="miRNA") 
  
  temp$stat <- temp %>% apply(1,function(x){
    temp.fish.cal <- fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE),alternative = "great")
    if (temp.fish.cal$estimate=="Inf") {
      temp.fish.cal.mod <- fisher.test(matrix(c(as.numeric(x['TSet']),as.numeric(x['BgSet'])-as.numeric(x['TSet']),as.numeric(x['TBg'])-as.numeric(x['TSet']),as.numeric(x['Bg'])-as.numeric(x['TBg'])-(as.numeric(x['BgSet'])-as.numeric(x['TSet']))),nrow = 2, byrow = TRUE)+1,alternative = "great")
      return(paste(temp.fish.cal.mod$estimate,temp.fish.cal$p.value,sep=":"))
    }else{
      return(paste(temp.fish.cal$estimate,temp.fish.cal$p.value,sep=":"))
    }
    
  })
  temp <- temp %>% separate(stat,c("od","pvalue"),sep = ":") %>% mutate(od=as.numeric(od),pvalue=as.numeric(pvalue)) %>% mutate(p_val_adj=p.adjust(pvalue,method="fdr"))
  return(temp)
}


FunPairedCor <- function(left_exp,right_exp,left_right,method='pearson') {
  
  left_right <- left_right %>% setNames(c("left","right"))
  temp <- left_exp %>% tibble::rownames_to_column("left") %>% tbl_df() %>% filter(left %in% left_right$left) %>% gather(sample,left_exp,-left) %>% inner_join(left_right,by="left",relationship = "many-to-many") %>% inner_join(right_exp %>% tibble::rownames_to_column("right") %>% tbl_df() %>% filter(right %in% left_right$right) %>% gather(sample,right_exp,-right),by=c("right","sample"),relationship = "many-to-many")
  
  temp.pos <- temp %>%  group_by(left,right)  %>% do(wtd_cor(.$left_exp,.$right_exp,weights = NULL,alternative="great",method=method)) 
  temp.neg <- temp %>%  group_by(left,right)  %>% do(wtd_cor(.$left_exp,.$right_exp,weights = NULL,alternative="less",method=method)) 
  temp.output <- temp.pos  %>% filter(r > 0) %>% bind_rows(temp.neg%>% filter(r < 0)  )
  return(temp.output)
}

FunPairedCor_noPvalue <- function(left_exp,right_exp,left_right,method='pearson') {
  
  left_right <- left_right %>% setNames(c("left","right"))
  temp <- left_exp %>% tibble::rownames_to_column("left") %>% tbl_df() %>% filter(left %in% left_right$left) %>% gather(sample,left_exp,-left) %>% inner_join(left_right,by="left",relationship = "many-to-many") %>% inner_join(right_exp %>% tibble::rownames_to_column("right") %>% tbl_df() %>% filter(right %in% left_right$right) %>% gather(sample,right_exp,-right),by=c("right","sample"),relationship = "many-to-many")
  
  temp.NP <-  temp %>%  group_by(left,right)  %>% do(wtd_cor(.$left_exp,.$right_exp,weights = NULL,method=method)) 
  rm(temp)
  gc()
  temp.output <- temp.NP %>% select(left,right,r) %>% mutate(r=ifelse(is.na(r),0,round(r,5)))
  return(temp.output)
}


Fun_ICM_TE_L1_regression <- function(sel_exp,meta,sel_embryo,tag) {
  temp.perf.list <- list()
  temp.coef.list <- list()
  
  for (em in sel_embryo ) {
    if (meta %>% filter(embryo==em) %>% pull(EML) %>% unique() %>% length() >1) {
      temp.test.cells <- meta %>% filter(embryo==em) %>% pull(cell)
      temp.train.cells <- meta$cell %>% setdiff(temp.test.cells)
      
      X.train <- sel_exp[,temp.train.cells] %>% t() 
      Y.train <- (meta %>% tibble::column_to_rownames("cell"))[temp.train.cells,"EML"] %>% factor(c("ICM","TE"))
      
      X.test <- sel_exp[,temp.test.cells]%>% t() 
      Y.test <- (meta %>% tibble::column_to_rownames("cell"))[temp.test.cells,"EML"]  %>% factor(c("ICM","TE")) 
      
      # L1 logistic regression
      set.seed(123)
      cvfit <- cv.glmnet(X.train, Y.train, family = "binomial", alpha = 1)#, nfolds = 10
      temp.model <- glmnet(X.train, Y.train, family = "binomial", alpha = 1, lambda = cvfit$lambda.min)
      
      # test
      temp.probs <- predict(temp.model, X.test, type = "response")[,1] %>% as.numeric()
      temp.pred <- predict(temp.model, X.test, type = "class")[,1] %>% as.vector()%>% factor(c("ICM","TE"))
      #print(cbind(temp.test.cells,temp.probs,as.vector(temp.pred),as.vector(Y.test)))
      
      # AUC
      temp.roc <- pROC::roc(Y.test,temp.probs)
      
      # cm calculation
      temp.cm <- caret::confusionMatrix(temp.pred,Y.test)
      
      #' output 
      temp.coef.list[[paste(em)]] <- coef(temp.model) %>% as.data.frame()%>% setNames("coef") %>% tibble::rownames_to_column("ID")  %>% tbl_df() %>% filter(ID!="(Intercept)") %>% mutate(tag=tag,sel_embryo=em) %>% filter(coef !=0)
      
      temp.perf.list[[paste(em)]] <- data.frame(sel_embryo=em,auc=as.numeric(temp.roc$auc),acc=as.numeric(temp.cm$overall["Accuracy"]),ba_acc=as.numeric(temp.cm$byClass["Balanced Accuracy"]),tag=tag)
    }
  }
    
  output <- list()
  output[["perf"]] <- temp.perf.list %>% do.call("bind_rows",.)
  output[["coef"]] <- temp.coef.list %>% do.call("bind_rows",.)
  return(output)
}
Fun_para_ICM_TE_L1_regression <- function(sel_exp_list,meta_list,perm_id,sel_embryo,st)  {
  sel_exp <- sel_exp_list[[st]]
  meta <- meta_list[[perm_id]]
  tag <-  paste(perm_id,st,sep=":")
  return(Fun_ICM_TE_L1_regression(sel_exp,meta,sel_embryo,tag)$perf)
}


Fun_para_miRNA_gene_cor <- function(perm_list,miRNA_exp,gene_exp,left_right,st,sel_cell=NULL)  {
  #tag <-  perm_list[[st]]
  left_data <- gene_exp
  left_data <- miRNA_exp
  
  
  colnames(left_data) <- perm_list[[st]][["cell_name"]]
  rownames(left_data) <- perm_list[[st]][["miRNA_name"]] 
  rownames(right_data) <- perm_list[[st]][["gene_name"]] 
  if (!is.null(sel_cell)) {
    left_data <- left_data[,sel_cell]
    right_data <- right_data[,sel_cell]
  }
  
  return(FunPairedCor_noPvalue(left_data ,right_data,left_right,method="spearman") %>% ungroup() %>% rename(miRNA=left,gene=right))
}

Fun_para_md_TF_cor <- function(perm_list,TF_exp,md_exp,left_right,st,sel_cell=NULL)  {
  #tag <-  perm_list[[st]]
  left_data <- TF_exp
  right_data <- md_exp
  
  colnames(left_data) <- perm_list[[st]][["cell_name"]]
  rownames(left_data) <- perm_list[[st]][["TF_name"]] 
  rownames(right_data) <- perm_list[[st]][["md_name"]] 
  if (!is.null(sel_cell)) {
    left_data <- left_data[,sel_cell]
    right_data <- right_data[,sel_cell]
  }
  
  return(FunPairedCor_noPvalue(left_data ,right_data,left_right,method="spearman") %>% ungroup() %>% rename(TF=left,md=right))
}