run_batch_analysis = function(data = data,compare_names = compare_names,org = "mmu"){
  
  if (!dir.exists(compare_names)) {
    dir.create(compare_names)
  }
  
  
  
  tempOutput = data
  names(tempOutput) = c("logFC","adj.P.Val")
  
  #------------------------------------------------------------------------------#
  # 筛选差异基因
  DEGs = tempOutput[which(tempOutput$adj.P.Val < 0.05 & tempOutput$logFC >= 2 | tempOutput$adj.P.Val < 0.05 & tempOutput$logFC <= -2),]
  AllDEGs = rownames(DEGs)
  length(AllDEGs)
  ## 火山图绘制
  
  data <-  tempOutput # 如果显著性不多 可以考虑用未校准的p值
  #data <- tempOutput[,c(1,4)] # 用未校准的p值
  # 绘制火山图
  library(ggplot2)
  library(ggrepel)
  
  colnames(data)
  
  names(data) <- c("log2FoldChange","padj")
  #data <- data[,c("log2FoldChange","padj")]
  data <- na.omit(data)
  #################
  # ggplot2绘制火山图
  data_tmp = data[which(data$log2FoldChange >= 2 | data$log2FoldChange <= -2),]
  data_tmp = data_tmp[order(data_tmp$padj),]

  #label = c(rownames(data_tmp)[1:10],rownames(data_tmp)[(nrow(data_tmp) - 9):nrow(data_tmp)])
  label = rownames(data_tmp)[1:20]
  
  
  #data <- data[order(data$log2FoldChange),]
  data$label = NA
  
  for (i in label) {
    data$label[which(rownames(data) == i)] = i
  }
  
  
  #data$label[which(rownames(data) %in% label)] = label
  #data$label[which(rownames(data) == "IL13")] = "IL13"
  #data$label[which(rownames(data) == "FGF19")] = "FGF19"
  
  data$regulate <- "None"
  data$regulate[which(data$log2FoldChange >= 2 & data$padj < 0.05)] = "UP"
  data$regulate[which(data$log2FoldChange <= -2 & data$padj < 0.05)] = "DOWN"
  data$regulate <- factor(data$regulate,levels = c("UP","None","DOWN"))
  
  table(data$regulate)
  #上调的有2448个，下调的有1103个
  #data$padj[which(data$padj == 0)] <- 1.196307e-302
  
  #data$padj[which(data$padj == 0)] = 2
  
  ylim = round(-log10(min(data$padj))) + 1
  
  data$padj[which(data$padj == 2)] = 10^(-ylim)
  
  ggplot(data,aes(log2FoldChange,-log10(padj),color = regulate)) +
    geom_point(alpha=0.4, size=2) +
    theme_bw(base_size = 12) +
    xlab("Log2(Fold change)") +
    ylab("-Log10(P.adj)") +
    theme(plot.title = element_text(size=15,hjust = 0.5)) +
    scale_colour_manual(values = c('red','gray',"blue")) +
    geom_hline(yintercept = -log10(0.05), lty = 4) +
    geom_vline(xintercept = c(-2, 2), lty = 4)+
    labs(title = "")+
    geom_label_repel(data = data, aes(label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",
                     show.legend = FALSE, max.overlaps = 10000) +
    ylim(c(0,ylim)) + ggtitle(label =  compare_names)
  table(data$regulate)
  
  
  
  ggsave(paste(compare_names,"/","Volcanoplot_",compare_names,".pdf",sep = ""),width = 8,height = 8)
  ggsave(paste(compare_names,"/","Volcanoplot_",compare_names,".png",sep = ""),width = 8,height = 8)
  write.csv(data, paste(compare_names,"/","DEG_sig_result",compare_names,".csv",sep = ""))
  
  library(scEasy)
  
  up = rownames(data[which(data$regulate == "UP"),])
  down = rownames(data[which(data$regulate == "DOWN"),])
  all = c(up,down)
  
  
  write(up,"./04_Relapsed/up.txt")
  write(down,"./04_Relapsed/down.txt")
  write(all,"./04_Relapsed/all.txt")
  
  dir.create(paste(compare_names,"/","Enrichment/","up",sep = ""),recursive = T)
  dir.create(paste(compare_names,"/","Enrichment/","down",sep = ""),recursive = T)
  dir.create(paste(compare_names,"/","Enrichment/","all",sep = ""),recursive = T)
  
  
  org_dir = getwd()
  setwd(paste(compare_names,"/","Enrichment/","up",sep = ""))
  scEasy::Run_Gene_Enrichment(geneList = up,org = org)
  
  setwd(org_dir)
  setwd(paste(compare_names,"/","Enrichment/","down",sep = ""))
  try(scEasy::Run_Gene_Enrichment(geneList = down,org = org))
  
  setwd(org_dir)
  setwd(paste(compare_names,"/","Enrichment/","all",sep = ""))
  scEasy::Run_Gene_Enrichment(geneList = all,org = org)
  
  setwd(org_dir)
}