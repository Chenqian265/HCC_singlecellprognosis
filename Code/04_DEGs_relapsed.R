library(Seurat)

sce = readRDS("./03_SCISSOR/01_scissor_result/sce_result_scissorCells.rds")

sce = subset(sce,scissor == "Scissor+")
table(sce$group)

Idents(sce) = "group"

dt = FindMarkers(sce,ident.1 = "Relapsed",ident.2 = "Primary",logfc.threshold = 0,min.pct = 0.1)
head(dt)

source("./Code/batch_analysis.R")

dt = dt[,c(2,1)]

names(dt) = c("logFC","adj.P.Val")

run_batch_analysis(data = dt,compare_names = "Scissor+ Relapsed vs Primary",org = "hsa")




DEGs = read.table("./03_SCISSOR/scissor+GENEs.txt")$V1

length(intersect(DEGs,all))


data = list(DEGs = all,markers = DEGs)

ggvenn::ggvenn(data = data,columns = c("DEGs","markers"))

ggsave("./04_Relapsed/venn.pdf",width = 8,height = 6)
ggsave("./04_Relapsed/venn.png",width = 8,height = 6)
