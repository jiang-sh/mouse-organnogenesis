#/software/biosoft/software/Rsoft/R3.6/bin/R

library(Seurat)
library(factoextra)
library(dendextend)
library(Matrix)

#x <- readRDS("/path/to/seurat/rds")
#Idents(x) <- x$sub_label
#cluster.average <- AverageExpression(x, return.seurat=TRUE)
#y <- t(cluster.average@assays$RNA@data)

y <- readRDS("subcluster.averages.0613.rds")

#Select highly variable genes
y <- FindVariableFeatures(y,selection.method="vst",nfeatures=10000)
varF <- y@assays$RNA@var.features
y <- t(y@assays$RNA@data[varF,])

#Compute distances and hierarchical clustering
dd <- get_dist(scale(y), method="euclidean")
hc <- hclust(dd, method="ward.D2")

#fviz_dend(hc,cex=0.5)
pdf("hc.378.top10000genes.euclidean.pdf",width=30,height=40)
fviz_dend(hc,k=60, cex=0.5, horiz=TRUE, color_labels_by_k = TRUE, ggtheme=theme_gray(),repel=TRUE)
dev.off()
