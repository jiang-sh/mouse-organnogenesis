library(ArchR)
packageVersion("ArchR")
addArchRThreads(threads = 12)
addArchRGenome("mm10")

path='/asnas/liujiang_group/yuhao/SingleCell_Mouse/01_Brain_cell'
setwd(path)
write(paste0(Sys.time()," START!"),file = "execute.log",append = T)

## 0 load and subset ##
load("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/Script_new/0_whole/TF_analysis/chromVAR/chromVAR_MM.RData")
brain_cell <- proj2@cellColData[which(as.character(proj2@cellColData$organ) %in% c("ForeBrain","MidBrain","HindBrain","SpinalCord")),]  ## 46770
brain <- proj2[rownames(brain_cell), ]
brain <- saveArchRProject(ArchRProj = brain, outputDirectory = "/asnas/liujiang_group/yuhao/SingleCell_Mouse/01_Brain_cell/ArchR_Brain", load = TRUE)
# saveArchRProject(ArchRProj = brain, outputDirectory = "/asnas/liujiang_group/yuhao/SingleCell_Mouse/01_Brain_cell/ArchR_Brain", load=TRUE,overwrite=TRUE)
write(paste0(Sys.time()," 0 Load success"),file = "execute.log",append = T)

## 1 Dimensionality Reduction with ArchR ##
brain <- addIterativeLSI(
  ArchRProj = brain,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)

write(paste0(Sys.time()," 1 Finish Dimensionality Reduction"),file = "execute.log",append = T)
# load("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/01_Brain_cell/ArchR_Brain/ArchR_dimensionality_reduction.RData")

## 2 Single-cell Embeddings ##
brain <- addUMAP(
  ArchRProj = brain,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)
write(paste0(Sys.time()," 2-1 Single-cell Embeddings"),file = "execute.log",append = T)

pdf("Plots/1_plotEmbedding_harmony.pdf")
p1 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "organ", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Clusters_sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "germ_celltype_ratio", embedding = "UMAP")
p1 + NoLegend()
p1
p2 + NoLegend()
p2
p3 + NoLegend()
p3
p4 + NoLegend()
p4

dev.off()
rm(list = c("p1","p2","p3", "p4"))
write(paste0(Sys.time()," 2-2 Single-cell Embeddings"),file = "execute.log",append = T)

## 2.1 Identify Marker Genes ##
markersGS <- getMarkerFeatures(
  ArchRProj = brain,
  useMatrix = "GeneScoreMatrix",
  groupBy = "germ_celltype_ratio",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList,'0_MarkerFeatures.csv')


## 3 label transfer ##
pdf("Plots/2_predictedScore_Co.pdf")
hist(
  brain@cellColData$predictedScore,
  xlab="prediction score",
  col="lightblue",
  xlim=c(0, 1),
  ylim=c(0, 2000),
  main="Kidney"
);
dev.off()
pdf("Plots/2_predictedScore_Co_UMAP.pdf")
plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "predictedScore", embedding = "UMAPHarmony")
#plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAPHarmony")
#plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAPHarmony")
dev.off()
write(paste0(Sys.time()," 3 Evaluate label transfer"),file = "execute.log",append = T)

save.image(file = "1_Save_Proj_Embedding.RData")

