rm(list = ls())
#---------------准备scissor输入数据------------------------------------------#
sc_dataset = readRDS("./02_Seurat_result/sce_final.rds")

metadata = readRDS("./00_Rawdata/prognostic/Rdata/survData_TCGA.rds")
metadata$samples = rownames(metadata)
metadata = metadata[,c(3,1,2)]
names(metadata)[2:3] = c("time","status")


bulk_dataset = readRDS("./00_Rawdata/prognostic/Rdata/exprSet_TCGA_tpm_tumor.rds")
bulk_dataset = as.matrix(bulk_dataset)

plot_title = "HCC"
name = "prognostic";alpha = 0.01;cutoff = 0.5



# devtools::install_github('sunduanchen/Scissor')
source("./Code/Run_scissor_pipline.R")
source("./Code/scissor_v1.R")
Run_scissor_pipline(sc_dataset = sc_dataset,metadata = metadata,bulk_dataset = bulk_dataset,
                    plot_title = plot_title,name = name,alpha = alpha,cutoff = cutoff)


# sc_dataset$scissor_ICR = "Not_defined"
# sc_dataset$scissor_ICR[which(sc_dataset$Condition == "Res+")] = "ICR_res"
# sc_dataset$scissor_ICR[which(sc_dataset$Condition == "Res-")] = "ICR_sen"
# UMAP_scissor <- DimPlot(sc_dataset, reduction = 'umap',
#                         group.by = 'scissor_ICR',
#                         cols = c('grey','royalblue','indianred1'),
#                         pt.size = 1, order = c("ICR_res","ICR_sen"))
# UMAP_scissor
# ggsave(paste(result_dir,"umap_scissor_prognosis_ICR.png",sep = "/"),width = 10,height = 10)
# ggsave(paste(result_dir,"umap_scissor_prognosis_ICR.pdf",sep = "/"),width = 10,height = 10)
# 
# 
# DimPlot(sc_dataset, reduction = 'umap',group.by = "scissor_ICR",split.by = 'celltype_main',ncol = 2,cols = c('indianred1','royalblue',"grey"),pt.size = 0.5)
# 
# ggsave(paste(result_dir,"umap_scissor_prognosis_ICR_splitby_samples.png",sep = "/"),width = 10,height = 10)
# ggsave(paste(result_dir,"umap_scissor_prognosis_ICR_splitby_samples.pdf",sep = "/"),width = 10,height = 10)

#saveRDS(sc_dataset,paste(result_dir,"sc_dataset_final_sceobj.rds",sep = "/"))

fs::dir_copy(path = "./Result/",new_path = "./03_SCISSOR/")


sce = readRDS("./03_SCISSOR/Result/sce_result_scissorCells.rds")
sc_dataset = sce
table(sce$scissor)

################################################################################
plot_scissor = function(sc_dataset = sc_dataset){
  ########################## KEGG and Go function analysis
  annotationDB = "org.Hs.eg.db"
  orgnm = "hsa"
  #Genelist <- read.csv(paste(result_dir,"DEG_result.csv",sep = "/"),row.names = 1)
  Idents(sc_dataset) <- "prognosis"
  Genelist <- FindMarkers(sc_dataset,ident.1 = "Unfavorable",ident.2 = "Favorable",logfc.threshold = 0,min.pct = 0.1,only.pos = TRUE)
  dt = Genelist
  
  
  Genelist <- Genelist[which(Genelist$avg_log2FC >= 2),]
  dim(Genelist)
  Genelist <- rownames(Genelist)
  library(scEasy)
  #source("./Run_Gene_Enrichment.R")
  Run_Gene_Enrichment(geneList = Genelist,org = "hsa")
  
  write(Genelist,"./03_SCISSOR/scissor+GENEs.txt")
  ########################################################
  dt <- FindAllMarkers(sc_dataset,logfc.threshold = 0,min.pct = 0.1,only.pos = TRUE)
  dt = dt[which(dt$cluster != "Not_defined"),]
  dt = dt[order(dt$cluster,dt$avg_log2FC,decreasing = T),]

  #saveRDS(dt,"DEGs_res_vs_sen.rds")
  write.csv(dt,"./03_SCISSOR/02_scissor+KEGG/DEGs_unfavirable_vs_favorable.csv")
  
  dt %>% group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% top_n(n = 5,wt = avg_log2FC ) -> dt_t
  
  sl = dt_t$gene
  
  
  
  
  library(scRNAtoolVis)
  jjVolcano(diffData = dt,base_size = 12,pSize = 2,flip = T,back.col = "lightblue",legend.position = c(0.95,0.9),log2FC.cutoff = 0.58,order.by = c("avg_log2FC"),myMarkers = sl)
  
 
  ggsave("./03_SCISSOR/02_scissor+KEGG/jjVolcano.pdf",width = 8.28,height = 6.52)
  ggsave("03_SCISSOR/02_scissor+KEGG/jjVolcano.png",width = 8.28,height = 6.52)
  # p5 = p4 + ggplot2::coord_flip()
  # p5
  # p5 <- p4 + ggplot2::scale_y_continuous(n.breaks = 6) +
  #       ggplot2::geom_label(ggplot2::aes(x = cluster,
  #                                        y = 0, label = cluster)) + ggplot2::theme(axis.line.y = ggplot2::element_blank(),
  #                                                                                  axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
  #       ggplot2::coord_flip()
  # }
  # p5
  
  setwd("./sl/")
  
  
  
  require(ggplot2)
  require(tidyverse)
  #require(enrichplot)
  gofile <- readxl::read_excel("./GO.xlsx")
  
  keggfile <- readxl::read_excel("./KEGG.xlsx")
  
  
  gofile <- as.data.frame(gofile)
  keggfile <- as.data.frame(keggfile)
  
  ##############################################################################
  go_res <- gofile
  # 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
  goBP <- subset(go_res,subset = (ONTOLOGY == "BP"))
  goCC <- subset(go_res,subset = (ONTOLOGY == "CC"))
  goMF <- subset(go_res,subset = (ONTOLOGY == "MF"))
  go.df <- rbind(goBP,goCC,goMF)
  
  go.df$ONTOLOGY <- factor(go.df$ONTOLOGY,levels = rev(c("BP","CC","MF")))
  go.df <- go.df[order(go.df$ONTOLOGY,go.df$Count,decreasing = T),]
  
  # 使画出的GO term的顺序与输入一致
  go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
  
  # 绘图
  go_bar <- ggplot(data = go.df, # 绘图使用的数据
                   aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
    geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
    coord_flip()+theme_classic()+ # 横纵坐标反转及去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
    #scale_fill_manual(values = rev(c("#26547C","#EF476F","#FFD166"))) +
    #scale_fill_manual(values  = c("#A13425","#E0C9A7","#283870")) +  # 配色方案 mumuxi推荐
    theme(axis.title = element_text(size = 13), # 坐标轴标题大小
          axis.text = element_text(size = 10), # 坐标轴标签大小
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
          legend.title = element_text(size = 13), # 图例标题大小
          legend.text = element_text(size = 11))# 图边距
  print(go_bar) + NoLegend()+theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))
  
  ggsave(go_bar,filename = "GO_Barplot_replot.pdf",width = 9,height = 7)
  ggsave(go_bar,filename = "GO_Barplot_replot.png",width = 9,height = 7)
  
  
  kegg_res <- keggfile
  
  
  kegg_res <- kegg_res[,colnames(kegg_res) != "ID"]
  #kegg_res <- kegg_res[order(kegg_res$qvalue,kegg_res$GeneRatio),]
  
  library(plyr)
  library(stringr)
  library(grid)
  library(ggplot2)
  e_data <- kegg_res
  #e_data <- e_data[,-1]
  # 分数转小数
  mixedToFloat <- function(x){
    x <- sapply(x, as.character)
    is.integer  <- grepl("^-?\\d+$", x)
    is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
    is.float <- grepl("^-?\\d+\\.\\d+$", x)
    is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
    stopifnot(all(is.integer | is.fraction | is.float | is.mixed))
    
    numbers <- strsplit(x, "[ /]")
    
    ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
           ifelse(is.float,    as.numeric(sapply(numbers, `[`, 1)),
                  ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                           as.numeric(sapply(numbers, `[`, 2)),
                         as.numeric(sapply(numbers, `[`, 1)) +
                           as.numeric(sapply(numbers, `[`, 2)) /
                           as.numeric(sapply(numbers, `[`, 3)))))
    
  }
  
  e_data_1 <- e_data
  e_data_1$GeneRatio = mixedToFloat(e_data_1$GeneRatio)
  e_data_1$Count = mixedToFloat(e_data_1$Count)
  log_name <- "-Log10(qvalue)"
  col_name_e_1 <- colnames(e_data_1)
  col_name_e_1 <- c(col_name_e_1,log_name)
  e_data_1$log_name <- log10(e_data_1$qvalue) * (-1)
  colnames(e_data_1) <- col_name_e_1
  
  e_data_1_freq <- as.data.frame(table(e_data_1$Description))
  colnames(e_data_1_freq) <- c("Description","ID")
  head(e_data_1_freq)
  
  
  e_data_2 <- merge(e_data_1,e_data_1_freq,by="Description")
  e_data_3 <- e_data_2[order(e_data_2$ID,
                             e_data_2$GeneRatio,
                             e_data_2$`-Log10(qvalue)`),]
  
  
  t_order <- unique(e_data_3$Description)
  e_data_1$Description <- factor(e_data_1$Description,
                                 levels = t_order,ordered = T)
  
  color_1 <- c("green","red")
  p <- ggplot(e_data_1,aes(x=GeneRatio,y=Description)) +
    labs(x="GeneRatio",y="GO description") + labs(title="")
  
  p
  p <- p + geom_point(aes(size=Count,color = `-Log10(qvalue)`)) +
    scale_color_gradient(low = color_1[1],high=color_1[2],name="-Log10(qvalue)")
  p
  p <- p + scale_y_discrete(labels=function(x) str_wrap(x,width = 60))
  p + theme_classic() + ggtitle("KEGG")+theme(plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
  
  ggsave(filename = "kegg_dotplot_replot.pdf",width = 9,height = 7)
  ggsave(filename = "kegg_dotplot_replot.png",width = 9,height = 7)
  print("The replot was finished!")
  
}
