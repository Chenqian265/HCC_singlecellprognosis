Run_scissor_pipline <- function(sc_dataset = sc_dataset,metadata = metadata,bulk_dataset = bulk_dataset,
                                plot = T,plot_title = "Epi_cluster",name = name,
                                alpha = alpha,cutoff = cutoff) {
  suppressPackageStartupMessages({
    library(Scissor)
    library(Seurat)
    library(tidyverse)
    library(org.Hs.eg.db)
  })
  col = CellChat::scPalette(20)
  if (!dir.exists("Result")) {
    dir.create("Result")
  }
  result_dir = "Result"
  
  UMAP_celltype = DimPlot(sc_dataset,reduction = "umap",label = T,pt.size = 1,label.size = 6,cols = col,rep = T) + theme_bw() + 
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  
  #UMAP_celltype
  
  if (identical(metadata$samples,colnames(bulk_dataset))) {
    if (colnames(metadata)[2] == "time") {
      print("Now run survival model")
      phenotype <- metadata[,2:3]
      colnames(phenotype) <- c("time", "status")
      
      infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha =  alpha,
                        family = "cox", cutoff =  cutoff,
                        Save_file = paste(result_dir,'Scissor_survival.RData',sep = "/"))
      # Reliability Significance Test
      load(paste(result_dir,'Scissor_survival.RData',sep = "/"))
      numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)#统计被Scissor选择的细胞数目
      print("SCISSOR选择的数目：")
      print(numbers)
      #result1 <- reliability.test(X, Y, network, alpha = 0.05, family = "cox", cell_num = numbers, n = 100, nfold = 10)
      #evaluate_summary <- evaluate.cell(paste(result_dir,'Scissor_survival.RData',sep = "/"), infos1, FDR_cutoff = 0.05, bootstrap_n = 100)
      
      # write out results
      #write.csv(evaluate_summary, paste(result_dir,'evaluate_summary.csv',sep = "/"))
      #write(paste("Reliability_test p value =",result1$p),paste(result_dir,'Reliability_test.txt',sep = "/"))
      
      # umap figures
      Scissor_select <- rep("none", ncol(sc_dataset))
      names(Scissor_select) <- colnames(sc_dataset)
      Scissor_select[infos1$Scissor_pos] <- "Scissor+"
      Scissor_select[infos1$Scissor_neg] <- "Scissor-"
      sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
      sc_dataset$scissor <- factor(sc_dataset$scissor,levels = c("Scissor+","Scissor-","none"))
      DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c("indianred1","royalblue","grey"), pt.size = 1)
      ggsave(paste(result_dir,"umap_scissor.png",sep = "/"))
      ggsave(paste(result_dir,"umap_scissor.pdf",sep = "/"))
      
      
      Scissor_select <- rep("none", ncol(sc_dataset))
      names(Scissor_select) <- colnames(sc_dataset)
      Scissor_select[infos1$Scissor_pos] <- "Unfavorable"
      Scissor_select[infos1$Scissor_neg] <- "Favorable"
      sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "prognosis")
      sc_dataset$prognosis <- factor(sc_dataset$prognosis,levels = c("Unfavorable","Favorable","none"))
      DimPlot(sc_dataset, reduction = 'umap', group.by = 'prognosis', cols = c("indianred1","royalblue","grey"), pt.size = 1)
      ggsave(paste(result_dir,"umap_scissor_prognosis.png",sep = "/"))
      ggsave(paste(result_dir,"umap_scissor_prognosis.pdf",sep = "/"))
      ## scissor cells
      sce_result <- subset(sc_dataset,scissor %in% c("Scissor+","Scissor-"))
      saveRDS(sce_result,paste(result_dir,"sce_result_scissorCells.rds",sep = "/"))
      saveRDS(sc_dataset,paste(result_dir,"sce_result_all.rds",sep = "/"))
      ### T test DEGs
      DefaultAssay(sce_result) <- "RNA"
      table(sce_result$scissor)
      Idents(sce_result) <- "scissor"
      DEG <- Seurat::FindMarkers(sce_result,ident.1 = "Scissor+",ident.2 = "Scissor-")
      DEG <- DEG[which(DEG$p_val_adj <= 0.05),]
      write.csv(DEG,paste(result_dir,"DEG_result.csv",sep = "/"))
    } else if (colnames(metadata)[2] == "status") {
      print("Now run  Classification Model")
      phenotype <- metadata$status
      names(phenotype) <- metadata$samples
      #identical(names(phenotype),colnames(bulk_dataset))
      tag <- c(0:(length(unique(metadata$status))-1))
      
      infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha =  alpha,
                        family = "binomial",
                        cutoff = cutoff,
                        Save_file = paste(result_dir,"Scissor_Classify_Model.RData",sep = "/"))
      
      
      # Reliability Significance Test
      load(paste(result_dir,'Scissor_Classify_Model.RData',sep = "/"))
      numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)#统计被Scissor选择的细胞数目
      #result1 <- reliability.test(X, Y, network, alpha = 0.05, family = "binomial", cell_num = numbers, n = 100, nfold = 10)
      #evaluate_summary <- evaluate.cell(paste(result_dir,'Scissor_Classify_Model.RData',sep = "/"), infos1, FDR_cutoff = 0.05, bootstrap_n = 100)
      
      #write out results
      #write.csv(evaluate_summary, paste(result_dir,'evaluate_summary.csv',sep = "/"))
      #write(paste("Reliability_test p value =",result1$p),paste(result_dir,'Reliability_test.txt',sep = "/"))
      
      # umap figures
      Scissor_select <- rep("none", ncol(sc_dataset))
      names(Scissor_select) <- colnames(sc_dataset)
      Scissor_select[infos1$Scissor_pos] <- "Scissor+"
      Scissor_select[infos1$Scissor_neg] <- "Scissor-"
      sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
      sc_dataset$scissor <- factor(sc_dataset$scissor,levels = c("Scissor+","Scissor-","none"))
      DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c("indianred1","royalblue","grey"), pt.size = 1)
      ggsave(paste(result_dir,"umap_scissor.png",sep = "/"),width = 10,height = 10)
      ggsave(paste(result_dir,"umap_scissor.pdf",sep = "/"),width = 10,height = 10)
      
      
      Scissor_select <- rep("none", ncol(sc_dataset))
      names(Scissor_select) <- colnames(sc_dataset)
      Scissor_select[infos1$Scissor_pos] <- paste(name,"+",sep = "")
      Scissor_select[infos1$Scissor_neg] <- paste(name,"-",sep = "")
      sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "Condition")
      
      sc_dataset$Condition <- factor(sc_dataset$Condition,levels = c(paste( name,"+",sep = ""),paste( name,"-",sep = ""),"none"))
      DimPlot(sc_dataset, reduction = 'umap', group.by = 'Condition', cols = c("indianred1","royalblue","grey"), pt.size = 1)
      ggsave(paste(result_dir,"umap_scissor_Condition.png",sep = "/"),width = 10,height = 10)
      ggsave(paste(result_dir,"umap_scissor_Condition.pdf",sep = "/"),width = 10,height = 10)
      ## scissor cells
      sce_result <- subset(sc_dataset,scissor %in% c("Scissor+","Scissor-"))
      saveRDS(sce_result,paste(result_dir,"sce_result_scissorCells.rds",sep = "/"))
      saveRDS(sc_dataset,paste(result_dir,"sce_result_all.rds",sep = "/"))
      ### T test DEGs
      DefaultAssay(sce_result) <- "RNA"
      table(sce_result$scissor)
      
      Idents(sc_dataset) <- "scissor"
      
      DEG <- Seurat::FindMarkers(sc_dataset,ident.1 = "Scissor+")
      DEG <- DEG[which(DEG$p_val_adj <= 0.05),]
      write.csv(DEG,paste(result_dir,"DEG_result.csv",sep = "/"))
    } else {
      stop("Please check you metadata colnames must be samples time status")
    }
  } else {
    stop("Please check you metadata and bulk data,Ensure consistent sample order!")
  }
  
  
  # 开始绘图
  print("Now plot the result")
  UMAP_scissor <- DimPlot(sc_dataset, reduction = 'umap',
                          group.by = 'scissor',
                          cols = c('grey','royalblue','indianred1'),
                          pt.size = 1, order = c("Scissor+","Scissor-"))
  UMAP_scissor
  ggsave(paste(result_dir,"umap_scissor_prognosis_new.png",sep = "/"),width = 10,height = 10)
  ggsave(paste(result_dir,"umap_scissor_prognosis_new.pdf",sep = "/"),width = 10,height = 10)
  
  
  patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scissor), ncol = 2)
  
  ggsave(paste(result_dir,"umap_scissor_prognosis_merge.png",sep = "/"),width = 14,height = 8)
  ggsave(paste(result_dir,"umap_scissor_prognosis_merge.pdf",sep = "/"),width = 14,height = 8)
  
  library(gplots)
  balloonplot(table(sc_dataset$scissor,sc_dataset$celltype))
  dev.copy2pdf(file = paste(result_dir,"balloonplot.pdf",sep = "/"),width = 14,height = 8)
}