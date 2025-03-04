meta = read.csv("./00_Rawdata/scdata/HCC_cell_metadata.csv",header = T)
rownames(meta) = meta$sample_name

exp = read.table("./00_Rawdata/scdata/exp_HCC_logTPM.txt",header = T)

library(tidyverse)
library(Seurat)
exp = column_to_rownames(exp,var = "gene")
exp = exp[,rownames(meta)]

sce = CreateSeuratObject(counts = exp,meta.data = meta)
sce@assays$RNA$data = sce@assays$RNA$counts


meta = sce@meta.data

celltype = str_split(meta$cell_type,pattern = "_",simplify = T)
celltype = celltype[,2]

meta$celltype = celltype

meta = meta[,c(1,2,3,7,8,9)]


sce@meta.data = meta


table(sce$HCC_type)
table(sce$tissue_source)
table(sce$orig.ident)


saveRDS(sce,"./00_Rawdata/scdata/sce_allcells.rds")


# 保留Tumor 样本

sce_T = subset(sce,tissue_source == "Tumor")
dim(sce_T)

table(sce_T$HCC_type)
unique(sce_T$orig.ident)
table(sce_T$orig.ident,sce_T$HCC_type)

saveRDS(sce_T,"./01_scRNAseq/sce.rds")
