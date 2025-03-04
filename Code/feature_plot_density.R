
plot_density <- function(obj,#单细胞seurat对象
                            marker,#需要展示的marker基因
                            dim=c("TSNE","UMAP"),
                            reduction = c(1,2), #1 = umap 2 = tsne
                            size,#点的大小
                            ncol=NULL#多个marker基因表达图排序
){
  suppressPackageStartupMessages({
    require(ggplot2)
    require(ggrastr)
    require(Seurat)
    require(patchwork)
  })
  
  #dim = "UMAP"
  #marker = c("CD3E","CCL5","DHX9")
  #obj = sce_all
  #size =1 
  
  #cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  cold <- colorRampPalette(c("#ede574","#008c9e","#005f6b"))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  if(dim=="TSNE"){
    
    xtitle = "tSNE1"
    ytitle = "tSNE2"
    
  }
  
  if(dim=="UMAP"){
    
    xtitle = "UMAP1"
    ytitle = "UMAP2"
  }
  
  
  if(length(marker)==1){
    reduction = c("umap","tsne")
    plot <- FeaturePlot(obj, features = marker,reduction = reduction[1])
    data <- plot$data
    
    
    if(dim=="TSNE"){
      
      colnames(data)<- c("x","y","ident","gene")
      
    }
    
    if(dim=="UMAP"){
      
      colnames(data)<- c("x","y","ident","gene")
    }
    
    
    
    p <- ggplot(data, aes(x, y)) +
      geom_point_rast(shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size) +
      geom_density_2d(data=data[data$gene>0,], 
                      aes(x=x, y=y), 
                      bins = 5, colour="black") +
      scale_fill_gradientn(colours = mypalette)+
      scale_colour_gradientn(colours = mypalette)+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0.5),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) 
    
    return(p)
    
  }else{
    
    gene_list <- list()
    
    
    
    for (i in 1:length(marker)) {
      reduction = c("umap","tsne")
      plot <- FeaturePlot(obj, features = marker[i],reduction = reduction[1])
      data <- plot$data
      
      
      if(dim=="TSNE"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      if(dim=="UMAP"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      gene_list[[i]] <- data
      names(gene_list) <- marker[i]
    }
    
    plot_list <- list()
    
    
    for (i in 1:length(marker)) {
      
      p <- ggplot(gene_list[[i]], aes(x, y)) +
        geom_point_rast(shape = 21, stroke=0.25,
                        aes(colour=gene, 
                            fill=gene), size = size) +
        geom_density_2d(data=gene_list[[i]][gene_list[[i]]$gene>0,], 
                        aes(x=x, y=y), 
                        bins = 5, colour="black") +
        scale_fill_gradientn(colours = mypalette)+
        scale_colour_gradientn(colours = mypalette)+
        theme_bw()+ggtitle(marker[i])+
        labs(x=xtitle, y=ytitle)+
        theme(
          plot.title = element_text(size=12, face="bold.italic", hjust = 0.5),
          axis.text=element_text(size=8, colour = "black"),
          axis.title=element_text(size=12),
          legend.text = element_text(size =10),
          legend.title=element_blank(),
          aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      plot_list[[i]] <- p
    }
    
    
    #Seurat::CombinePlots(plot_list, ncol = ncol)
    
    wrap_plots(plot_list, ncol = ncol)
  }
  
  
}
