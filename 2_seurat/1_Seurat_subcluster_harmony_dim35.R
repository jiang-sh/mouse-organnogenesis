library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(biomaRt)
library(umap)

pca = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/subcluster/PCA/Seurat_wholeMap_harmony_"
mark = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/subcluster/test/marker/Seurat_wholeMap_harmony"
rds = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/subcluster/test/rds/Seurat_wholeMap_harmony"
pdf.file = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/subcluster/test/pdf/Seurat_wholeMap_harmony"


#all <- readRDS("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/wholeMap/Seurat_wholeMap_harmony_cell_runharmony_r0.6.rds")

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- read.table("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/cellcycle_s.txt")
g2m.genes <- read.table("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/cellcycle_g2m.txt")

############  subcluster each main cluster #################################
dimension = 35
cluster.ls = list()
for (i in c(7,32)){
  for (j in c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1)){
    name = paste0("cluster", i)
    cluster.ls[[name]] <- readRDS(paste0(pca, name, "_PCA.rds"))
    
    ############### 7_Cluster the cells #############################
    # use appropriate dim and find appropriate resolution
    #cluster.ls[[name]] <- FindNeighbors(cluster.ls[[name]], dims = 1:dimension)
    #cluster.ls[[name]] <- FindClusters(cluster.ls[[name]], resolution = j, algorithm = 1)
    
    ############# 8_Run non-linear dimensional reduction (UMAP/tSNE) ################
    #cluster.ls[[name]] <- RunUMAP(cluster.ls[[name]], dims = 1:dimension)

    ############### 7_Run harmony ###################################
    options(repr.plot.height = 2.5, repr.plot.width = 6)
    pdf(paste0(pdf.file, "_", name,  "_dim", dimension, "r", j, "_05_harmony_reduction.pdf"))
    cluster.ls[[name]] <- cluster.ls[[name]] %>%
      RunHarmony("batch", plot_convergence = TRUE)

    ### harmony embeddings ###
    harmony_embeddings <- Embeddings(cluster.ls[[name]], 'harmony')
    harmony_embeddings[1:5, 1:5]

    print(DimPlot(object = cluster.ls[[name]], reduction = "harmony", pt.size = .1, group.by = "batch", do.return = TRUE))
    print(DimPlot(cluster.ls[[name]], reduction = "harmony", group.by = "sample", pt.size = .1))
    print(DimPlot(cluster.ls[[name]], reduction = "harmony", group.by = "germlayer", pt.size = .1))
    dev.off()

    ########## visualization ########################### 
    cluster.ls[[name]] <- cluster.ls[[name]] %>%
      RunUMAP(reduction = "harmony", dims = 1:dimension) %>%
      FindNeighbors(reduction = "harmony", dims = 1:dimension) %>%
      FindClusters(resolution = j) %>%
      identity()
    

    pdf(paste0(pdf.file, "_", name, "_dim", dimension, "r", j, "_05_harmony_cluster.pdf"))
    options(repr.plot.height = 4, repr.plot.width = 6)
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = TRUE, pt.size = .1))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = T, group.by = "sample"))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = F, group.by = "sample"))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = T, group.by = "germlayer"))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = T, group.by = "Phase"))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = T, group.by = "batch"))
    print(DimPlot(cluster.ls[[name]], reduction = "umap", label = T, group.by = "DF_results"))
    dev.off()
    
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    name.markers <- FindAllMarkers(cluster.ls[[name]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(name.markers, file = paste0(mark, "_", name, "_dim", dimension, "r", j, "_clusters.markers.xls"), sep = "\t", quote = F, col.names = T, row.names = T)
    saveRDS(cluster.ls[[name]], file = paste0(rds, "_", name, "_dim", dimension, "r", j, "_cluster.rds")) 
  }
}
