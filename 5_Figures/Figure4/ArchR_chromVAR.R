library(ArchR)
library(Seurat)
library(Matrix)

addArchRThreads(threads = 8)
addArchRGenome("mm10")

setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/Script_new/0_whole/TF_analysis/chromVAR")
proj2 <- loadArchRProject("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/whole_Tss12000_v3/Save-Proj5_chromVAR_2")
### 2. chromVAR ###
if("Motif" %ni% names(proj2@peakAnnotation)){
    proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
}

proj2 <- addBgdPeaks(proj2)

### 2.1  compute per-cell deviations ###
proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "Motif",
  force = TRUE
)

MM <- getMatrixFromProject(
  ArchRProj = proj2,
  useMatrix = "MotifMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
saveRDS(MM,file = "chromVAR_Motif_matirx.rds")

MM@assays@data@listData$deviations
MM@assays@data@listData$z

### 2.2 getVarDeviations ###
plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)
pdf("0_plot_motif_Deviations.pdf")
plot(plotVarDev)
dev.off()

### select part of motif to plot ###
motifs <- c("Fox2a", "Pou3f3")
markerMotifs <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
### grep z score ###
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% c("z:Foxa2_307", "z:Pou3f3_540")]
markerMotifs

## expression pattern ##
markerRNA <- c( "Fox2a", "Pou3f3")

p2 <- plotEmbedding(
    ArchRProj = proj2,
    colorBy = "GeneIntegrationMatrix",
    name = sort(markerRNA),
    embedding = "UMAPHarmony",
    continuousSet = "blueYellow",
    imputeWeights = getImputeWeights(proj2)
)

p3 <- plotEmbedding(
    ArchRProj = proj2,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(proj2)
)


for (i in 1:2){
  pdf(paste0(as.character(markerRNA[i]),"_expression_chromVAR_embedding.pdf"))
  print(p2[[as.character(markerRNA[i])]])
  print(p3[[as.character(markerMotifs[i])]])
  dev.off()
}


