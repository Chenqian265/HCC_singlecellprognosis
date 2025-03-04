library(Seurat)
library(harmony)

sce = readRDS("./01_scRNAseq/sce.rds")
DefaultAssay(sce) = "RNA"
sce <- FindVariableFeatures(sce,selection.method = "vst",nfeatures = 2000)
sce <- ScaleData(sce)
sce <- RunPCA(sce,features = VariableFeatures(sce),reduction.name = "pca")
sce <- RunHarmony(sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
sce <- RunUMAP(sce,reduction = "harmony",dims = 1:20,reduction.name = "umap")

ElbowPlot(sce,reduction = "pca")
export::graph2pdf(file = "./02_Seurat_result/ElbowPlot_PCA.pdf",width = 8,height = 8)
export::graph2png(file = "./02_Seurat_result/ElbowPlot_PCA.png",width = 8,height = 8)

ElbowPlot(sce,reduction = "harmony")
export::graph2pdf(file = "./02_Seurat_result/ElbowPlot_harmony.pdf",width = 8,height = 8)
export::graph2png(file = "./02_Seurat_result/ElbowPlot_harmony.png",width = 8,height = 8)


sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20)
### 设置多个resolution选择合适的resolution
sce <- FindClusters(sce, resolution = 0.5)

Idents(sce) <- "celltype"
sce$celltype[which(sce$celltype %in% c("Mye.","Mey."))] = "Mye"


sce$group = sce$HCC_type


DimPlot(sce,reduction = "umap",label = T)


source("./Code/source.R")




library(ggplot2)
library(Seurat)
library(presto)
library(scRNAtoolVis)

dest_dir = "./02_Seurat_result/"

col = Col_list(length(unique(sce$celltype)))
DimPlot(sce,reduction = "umap",label = T,pt.size = 0.8,label.size = 8,cols = col,repel = T) + theme_bw() + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))
ggsave(paste0(dest_dir,"/umap_cluster_all.pdf"),width = 10,height = 10)
ggsave(paste0(dest_dir,"/umap_cluster_all.png"),width = 10,height = 10)
dev.off()


DimPlot(sce,reduction = "umap",label = T,pt.size = 0.8,label.size = 6,cols = col,split.by = "group",repel = T) + theme_bw() + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))
ggsave(paste0(dest_dir,"/umap_cluster_splitby_sample.pdf"),width = 10,height = 5)
ggsave(paste0(dest_dir,"/umap_cluster_splitby_sample.png"),width = 10,height = 5)
dev.off()


DefaultAssay(sce) <- "RNA"
df = FindAllMarkers(sce,logfc.threshold = 0,only.pos = F,min.pct = 0.1)

jjVolcano(diffData = df,aesCol = c('purple','orange'),tile.col = col,size = 3,base_size = 18,
          back.col = "lightgrey",topGeneN = 3,pSize = 1,legend.position = c(0.8,0.95))

ggsave(paste0(dest_dir,"/volcanoplot_cluster_markers.pdf"),width = 14,height = 5)
ggsave(paste0(dest_dir,"/volcanoplot_cluster_markers.png"),width = 14,height = 5)
dev.off()

# plot
df$ptFC = df$pct.1 - df$pct.2

markers = df %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(ptFC > 0.5) %>% slice_head(n = 5)

Seurat::DoHeatmap(subset(sce,downsample = 200),markers$gene)
export::graph2pdf(file = paste0(dest_dir,"/Top_Doheatmap.pdf"),width = 14,height = 14)
export::graph2png(file = paste0(dest_dir,"/Top_Doheatmap.png"),width = 14,height = 14)

write.csv(df,file = paste0(dest_dir,"/allmarkers_cluster.csv"),row.names = T)
write.csv(markers,file = paste0(dest_dir,"/Top5_markers_cluster.csv"),row.names = T)





col =  rev(col)


dat <- table(sce$celltype,sce$group)
dat <- prop.table(dat,margin = 1)
dat
write.csv(dat,"./Ratio_per_celltype.csv")


dat <- read.csv("./Ratio_per_celltype.csv")
names(dat)[1] = "celltype"
dat


x_levels <- unique(sce$celltype)
fill.colors <- col
width<-  0.7
dat <- dat %>% pivot_longer(!celltype)
names(dat)[2] = "group"

dat%>%
  mutate(group=factor(group)) %>%
  mutate(celltype=factor(celltype,levels = rev(x_levels))) %>%
  ggplot(aes(x=value,y=group))+
  geom_bar(aes(fill= celltype),
           stat="identity",
           position = "fill",
           width = width)+
  scale_fill_manual(values = col) + theme_classic() -> p1

p1

p1+     theme_classic()+
  scale_x_continuous(position = "top",
                     expand = expansion(mult=c(0,0)),
                     breaks = seq(0,1,by=0.1),
                     labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0))+
  theme(axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(face="bold",size=10),
        legend.title = element_blank())+
  labs(title="Percentage of cells per celltype")
ggsave("Ratio_per_celltype.pdf",width = 10,height = 4)
ggsave("Ratio_per_celltype.png",width = 10,height = 4)

#------------------------------------------------------------------------------#
dat <- table(sce$group,sce$celltype)
dat <- prop.table(dat,margin = 2)
dat
write.csv(dat,"./Ratio_per_group.csv")


dat <- read.csv("./Ratio_per_group.csv")
names(dat)[1] = "group"
dat


x_levels <- unique(sce$group)
fill.colors <- c("#ff7473","#ffc952","#47b8e0","#84bd00")


width<-  0.5

dat %>%
  pivot_longer(!group) %>%
  mutate(name=factor(name)) %>%
  mutate(group=factor(group,levels = rev(x_levels))) %>%
  ggplot(aes(x = value,y = name))+
  geom_bar(aes(fill=group),
           stat="identity",
           position = "fill",
           width = width)+
  scale_fill_manual(values = rev(fill.colors)) + theme_classic() -> p1

p1

p1 + theme_classic()+
  scale_x_continuous(position = "top",
                     expand = expansion(mult=c(0,0)),
                     breaks = seq(0,1,by=0.1),
                     labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0))+
  theme(axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(),
        axis.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(face="bold",size=10),
        legend.title = element_blank())+
  labs(title="Percentage of cells per group")
ggsave("Ratio_per_group.pdf",width = 10,height = 4)
ggsave("Ratio_per_group.png",width = 10,height = 4)




sce <- ScaleData(sce,features = rownames(sce))


Idents(sce) <- "celltype"
Allmarkers <- df

Allmarkers = transform(Allmarkers,FC.pt = (pct.1 - pct.2))


#dir.create("./Result")
write.csv(Allmarkers,"./02_Seurat_result//Allmarkers_anno.csv")


Allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(FC.pt > 0.6) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

write.csv(top10,"./02_Seurat_result//Allmarkers_anno_top10.csv")


Allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(FC.pt > 0.6) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

write.csv(top5,"./02_Seurat_result//Allmarkers_anno_top5.csv")





pdf("./02_Seurat_result//Pheatmap_top10.pdf",height = 20,width = 16)
p = DoHeatmap(sce, features = top10$gene) + NoLegend()
print(p)
dev.off()


pdf("./02_Seurat_result//Pheatmap_top5.pdf",height = 20,width = 16)
p = DoHeatmap(sce, features = top5$gene) + NoLegend()
print(p)
dev.off()




png("./02_Seurat_result//Pheatmap_top10.png",height = 2000,width = 1600)
p = DoHeatmap(sce, features = top10$gene) + NoLegend()
print(p)
dev.off()


png("./02_Seurat_result//Pheatmap_top5.png",height = 2000,width = 1600)
p = DoHeatmap(sce, features = top5$gene) + NoLegend()
print(p)
dev.off()


#############################################################################
source("./Code/feature_plot_density.R")
genelist <- top5$gene
if (!dir.exists("./02_Seurat_result//Featureplot/")) {
  dir.create("./02_Seurat_result//Featureplot")
}
for(i in genelist) {
  #i = "Gas6"
  #source("../code/feature_plot_density.R")
  plot_density(sce,i,dim = "UMAP",reduction = 1,size = 1)
  ggsave(paste("./02_Seurat_result//Featureplot//plotdensity_",i,".pdf",sep = ""),width = 10,height = 10)
  ggsave(paste("./02_Seurat_result//Featureplot//plotdensity_",i,".png",sep = ""),width = 10,height = 10)
}

Allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3


gl <- unique(top3$gene)

sce_sub <- sce[gl,]

sl = gl[gl %in% rownames(sce_sub)]


p = DotPlot(sce_sub, features = sl) +
  theme_bw()+ theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 0.2,vjust=0.5,angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
print(p)
ggsave("./02_Seurat_result//allcells_Dotplot_top3.pdf",width = 14,height = 10)
ggsave("./02_Seurat_result//allcells_Dotplot_top3.png",width = 14,height = 10)



#-----------------------------------------------------------------------------#
gl = c("CD3E","CD19","AKR1C4","CD14")
gl


if (!dir.exists("02_Seurat_result/KeyGenes/")) {
  dir.create("02_Seurat_result/KeyGenes/")
}
for(i in gl) {
  #i = "Gas6"
  #source("../code/feature_plot_density.R")
  plot_density(sce,i,dim = "UMAP",reduction = 1,size = 1)
  ggsave(paste("./02_Seurat_result/KeyGenes//plotdensity_",i,".pdf",sep = ""),width = 10,height = 10)
  ggsave(paste("./02_Seurat_result/KeyGenes///plotdensity_",i,".png",sep = ""),width = 10,height = 10)
}




gl = c("CD3D","NKG7","CD79A","JCHAIN","LYZ","FCER1G","CD14","DCN","AKR1C4")

plot_density(sce,gl,dim = "UMAP",reduction = 1,size = 1)
ggsave("./02_Seurat_result//Mainfeatureplot.pdf",width = 12,height = 12)
ggsave("./02_Seurat_result//Mainfeatureplot.png",width = 12,height = 12)

#-----------------------------------------------------------------------------#
# gl = c("RB1","VCP","MYCN","ZEB1","CXCL13")
# gl
# 
# source("./Code/source.R")
# sce$group = factor(sce$group,levels = c("Normal","Tumor"))
# VlnPlot(sce,features = gl,split.by = "group",cols = Col_list(2))
# 
# ggsave("./18_scRNAseq/04_KeyGenes//Mainvlnplot.pdf",width = 12,height = 8)
# ggsave("./18_scRNAseq/04_KeyGenes//Mainvlnplot.png",width = 12,height = 8)
# 
# 
# exp = AggregateExpression(sce,features = gl,group.by = c("celltype","group"),assays = "RNA")
# exp = exp$RNA
# exp = as.matrix(exp)
# 
# pheatmap::pheatmap(exp,scale = "row",cluster_cols = F)
# export::graph2pdf(file = "./18_scRNAseq/04_KeyGenes/pheatmap_key.pdf",width = 12,height = 9)
# 
# saveRDS(sce,"./18_scRNAseq/sce_final.rds")



col = rev(col)
DimPlot(sce,reduction = "umap",label = T,pt.size = 0.8,label.size = 4,cols = col,split.by = "orig.ident",repel = T,ncol = 3) + theme_bw() + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))


ggsave(paste0(dest_dir,"/umap_cluster_splitby_orig.ident.pdf"),width = 11,height = 12)
ggsave(paste0(dest_dir,"/umap_cluster_splitby_orig.ident.png"),width = 11,height = 12)
dev.off()

saveRDS(sce,"./02_Seurat_result/sce_final.rds")








