source("./Code/feature_plot_density.R")
library(Seurat)
gl = read.table("./05_Prognostic/cox/07.multiCox.xls",header = T)
gl = gl$id

sce = readRDS("./02_Seurat_result/sce_final.rds")

plot_density(sce,gl,dim = "UMAP",reduction = 1,size = 1)
ggsave("./05_Prognostic/modelgene_featureplot.pdf",width = 9.5,height = 6.5)
ggsave("./05_Prognostic/modelgene_featureplot.png",width = 9.5,height = 6.5)



