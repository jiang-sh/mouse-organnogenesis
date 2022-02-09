#Xqlogin -q core40 -l nodes=1:ppn=12,mem=90gb,walltime=900:00:00

library('ArchR')
library('dplyr')
library('ComplexHeatmap')
setwd('/xtdisk/jiangl_group/yuchw/projects/Mouse_Organogenesis_v3/Gonad_Development/Gonad_Development_v2/Gonad_ArchR/')

addArchRThreads(threads = 12)
addArchRGenome("mm10")

proj_sc <- loadArchRProject("SupportingCell_ArchRProject")

fragments.files=c()
sample.names=c("Gonad-male_E11_5","Gonad-male_E12_5","Gonad-male_E13_5",
	"Gonad-female_E11_5","Gonad-female_E12_5","Gonad-female_E13_5");

for (index in 1:length(sample.names)){
  name_i = sample.names[index]
  fragments.files[index]=gsub("organ",name_i,"/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_barcode_multiplet_ArchR_bap2_v2/scATAC_organ_bap/final/scATAC_organ.fragments.tsv.gz")
}

## 1_Create Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles=fragments.files,
  sampleNames=paste0("scATAC_", sample.names),
  filterTSS=0,
  filterFrags=0,
  addTileMat=TRUE,
  addGeneScoreMat=TRUE
)

ArrowFiles

## 2_Doublet Finding
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

## 3.1_Creating An ArchRProject
proj <- ArchRProject(
  ArrowFiles=ArrowFiles,
  outputDirectory="Gonad_ArchRProject",
  copyArrows=TRUE
)

cellid <- read.table("Gonad_QC_rmD_barcode.txt", comment.char="")
proj <- proj[cellid$V1,]

proj$sample <- proj$Sample
proj$sample <- gsub("scATAC_Gonad-female_E11_5","E11-GF",proj$sample)
proj$sample <- gsub("scATAC_Gonad-female_E12_5","E12-GF",proj$sample)
proj$sample <- gsub("scATAC_Gonad-female_E13_5","E13-GF",proj$sample)
proj$sample <- gsub("scATAC_Gonad-male_E11_5","E11-GM",proj$sample)
proj$sample <- gsub("scATAC_Gonad-male_E12_5","E12-GM",proj$sample)
proj$sample <- gsub("scATAC_Gonad-male_E13_5","E13-GM",proj$sample)

proj$gender <- proj$Sample
proj$gender <- gsub("scATAC_Gonad-female_E11_5","female",proj$gender)
proj$gender <- gsub("scATAC_Gonad-female_E12_5","female",proj$gender)
proj$gender <- gsub("scATAC_Gonad-female_E13_5","female",proj$gender)
proj$gender <- gsub("scATAC_Gonad-male_E13_5","male",proj$gender)
proj$gender <- gsub("scATAC_Gonad-male_E12_5","male",proj$gender)
proj$gender <- gsub("scATAC_Gonad-male_E11_5","male",proj$gender)

proj$stage <- proj$Sample
proj$stage <- gsub("scATAC_Gonad-female_E11_5","E11.5",proj$stage)
proj$stage <- gsub("scATAC_Gonad-female_E12_5","E12.5",proj$stage)
proj$stage <- gsub("scATAC_Gonad-female_E13_5","E13.5",proj$stage)
proj$stage <- gsub("scATAC_Gonad-male_E13_5","E13.5",proj$stage)
proj$stage <- gsub("scATAC_Gonad-male_E12_5","E12.5",proj$stage)
proj$stage <- gsub("scATAC_Gonad-male_E11_5","E11.5",proj$stage)

table(proj$sample)
#E11-GF E11-GM E12-GF E12-GM E13-GF E13-GM 
#  7228   6615  10538  10557  10547   9460 
#  7139   6568  10418  10406  10525   9381

## 3.2_Filter Doublets
proj <- filterDoublets(proj)
table(proj$sample)
#E11-GF E11-GM E12-GF E12-GM E13-GF E13-GM 
#  6706   6178   9428   9552  10075   8566 

## 3.3_Plotting Sample Statistics from an ArchRProject
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj,
    groupBy = "sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics-whole.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)


## 4_Dimensionality Reduction with ArchR
proj2 <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed = 2,
    sampleCellsPre=30000,
    force = TRUE
)
# Run Harmony - batch effect correlation
proj2 <- addHarmony(
    ArchRProj=proj2,
    reducedDims="IterativeLSI",
    name="Harmony",
    groupBy="stage",
    force=TRUE
)
#Clustering
proj2 <- addClusters(
    input=proj2,
    reduceDims="Harmony",
    method="Seurat",
    name="Clusters",
    resolution=0.8
)

saveArchRProject(
    ArchRProj=proj2,
    outputDirectory="Gonad_ArchRProject-beforeUMAP",
    overwrite=TRUE,
    load=T,
    dropCells = FALSE,
    logFile=createLogFile("saveArchRProject"),
    threads=10
)

proj2 <- addUMAP(
    ArchRProj = proj,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 20,
    minDist = 0.6,
    metric = "cosine",
    force=TRUE
)

saveArchRProject(
    ArchRProj=proj2,
    outputDirectory="Gonad_ArchRProject-afterUMAP",
    overwrite=TRUE,
    load=T,
    dropCells = FALSE,
    logFile=createLogFile("saveArchRProject"),
    threads=12
)

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4,name="Plot-UMAP2Harmony-Sample-Clusters.0128.pdf",ArchRProj = proj2, addDOC = FALSE, width = 10, height = 5)

## Identifying Marker Features
markersGS <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerGenes <- c(
    "Col1a1", "Col1a2", "Ptn",       #Stromal Progenitor Cell
    "Insl3", "Cyp17a1", "Cyp11a1",   #Fetal Leydig Cell
    "Upk3b", "Myl7", "Krt7",         #Mesothelial Cell
    "Lgals7", "Kctd14", "Fst", "Wnt4", #Pre-Granulosa Cell
    "Ly6a", "Clec18a", "Lrp2",       #bi-potential state
    "Amh", "Cst9", "Ptgds",          #Sertoli Cell
    "Dppa5a", "Dppa3", "Pou5f1",     #Primordial Germ Cell
    "Sycp3", "Tex101", "Tex12",      #Meotic Germ Cell
    "Egfl7", "Cd34", "Plvap",        #Endothelial Cell
    "C1qb", "C1qc", "Cxcl2", 	     #Macrophage
    "Hba-x", "Hbb-y", "Hbb-bs"       #Erythroblast 
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot",annotation_legend_side="bot")
plotPDF(heatmapGS, name="GeneScores-Marker-Heatmap",width=10,height=8,ArchRProj=proj2, addDOC=FALSE)

## Integration
# Unconstrained Integration
proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI2",
    seRNA = gonad,
    addToArrow = FALSE,
    groupRNA = "celltype2",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# Constrained Integration
cM <- as.matrix(confusionMatrix(proj2$Clusters, proj2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1, which.max)]
cbind(preClust,rownames(cM))

cGerm <- "08|09|10|11"
cSupporting <- "01|02|22|23|24"
cStromal <- "05|06|07|13|14|15|17|18|19|20|21"
cCommon <- "03|04|12|16"

clustGerm <- c("C8","C9","C10","C11")
clustSupporting <- c("C1","C2","C22","C23","C24")
clustStromal <- c("C5","C6","C7","C13","C14","C15","C17","C18","C19","C20","C21")
clustCommon <- c("C3","C4","C12","C16")

rnaGerm <- rownames(gonad@meta.data)[grep("Germ_Cell",gonad$celltype2)]
rnaSupporting <- rownames(gonad@meta.data)[grep("Granulosa|Sertoli|Supporting",gonad$celltype2)]
rnaStromal <- rownames(gonad@meta.data)[grep("Stromal|Mesothelial|Leydig",gonad$celltype2)]
rnaCommon <- rownames(gonad@meta.data)[grep("Erythroblast|Macrophage|Endothelial",gonad$celltype2)]

groupList <- SimpleList(
    Germ = SimpleList(
        ATAC = proj3$cellNames[proj3$Clusters5 %in% clustGerm],
        RNA = rnaGerm
    ),
    Supporting = SimpleList(
       	ATAC = proj3$cellNames[proj3$Clusters5 %in% clustSupporting],
       	RNA = rnaSupporting
    ),
    Stromal = SimpleList(
       	ATAC = proj3$cellNames[proj3$Clusters5 %in% clustStromal],
       	RNA = rnaStromal
    ),
    Common = SimpleList(
       	ATAC = proj3$cellNames[proj3$Clusters5 %in% clustCommon],
       	RNA = rnaCommon
    )
)

proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI2",
    seRNA = gonad,
    addToArrow = TRUE, 
    groupList = groupList,
    groupRNA = "celltype2",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

pal <- paletteDiscrete(values = gonad$celltype2)

p1 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Un",
    pal = pal,
    embedding = "UMAPHarmony"
)
p2 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedScore_Un",
    embedding = "UMAPHarmony"
)
plotPDF(p1, p2, name = "Plot-UMAP-Un-Integration.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

p3 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Co",
    pal = pal,
    embedding = "UMAPHarmony"
)
p4 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedScore_Co",
    embedding = "UMAPHarmony"
)
plotPDF(p3, p4, name = "Plot-UMAP-Co-Integration.2.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

#plot prediction score  histogram
pdf("Gonad_ArchRProject-afterLT/Plots/LT_Prediction_Score.0131.pdf", width=8, height=5)
hist(
    proj3$predictedScore_Un,
    xlab="prediction score",
    col="lightblue",
    xlim=c(0, 1),
    main="scATAC unconstrained integration prediction score"
);
abline(v=0.4, col="red", lwd=2, lty=2);
table(proj3$predictedScore_Un > 0.4);

hist(
    proj3$predictedScore_Co,
    xlab="prediction score",
    col="lightblue",
    xlim=c(0, 1),
    main="scATAC constrained integration prediction score"
);
abline(v=0.4, col="red", lwd=2, lty=2);
table(proj3$predictedScore_Co > 0.4);
dev.off()

proj4 <- addGroupCoverages(ArchRProj=proj3, groupBy="predictedGroup_Co")
pathToMacs2 <- findMacs2()
proj4 <- addReproduciblePeakSet(
    ArchRProj = proj4, 
    groupBy = "predictedGroup_Co", 
    pathToMacs2 = pathToMacs2
)

getPeakSet(proj4)

proj5 <- addPeakMatrix(proj4)
getAvailableMatrix(proj5)

saveArchRProject(
    ArchRProj=proj5,
    outputDirectory="Gonad_ArchRProject-afterPeakCalling",
    overwrite=TRUE,
    load=T,
    dropCells = FALSE,
    logFile=createLogFile("saveArchRProject"),
    threads=12
)

####### second round LSI with PeakMatrix #######
####### third round LSI with GeneScoreMatrix #######

## 4_Dimensionality Reduction with ArchR
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed = 2,
    sampleCellsPre=40000,
    force = TRUE
)
# Run Harmony - batch effect correlation
proj <- addHarmony(
    ArchRProj=proj,
    reducedDims="IterativeLSI",
    name="Harmony",
    groupBy="stage",
    force=TRUE
)
#Clustering
proj3 <- addClusters(
    input=proj3,
    reducedDims="Harmony2",
    method="Seurat",
    name="Clusters2",
    resolution=0.8,
    force=TRUE
)
proj3 <- addClusters(
    input=proj3,
    reducedDims="Harmony",
    method="Seurat",
    name="Clusters",
    resolution=0.8,
    force=TRUE
)
proj3 <- addClusters(
    input=proj3,
    reducedDims="Harmony3",
    method="Seurat",
    name="Clusters3",
    resolution=0.8,
    force=TRUE
)
proj3 <- addClusters(
    input=proj3,
    reducedDims="",
    method="Seurat",
    name="Clusters4",
    resolution=0.8,
    force=TRUE
)
proj3 <- addClusters(
    input=proj3,
    reducedDims="IterativeLSI2",
    method="Seurat",
    name="Clusters5",
    resolution=0.8,
    force=TRUE
)
proj3 <- addClusters(
    input=proj3,
    reducedDims="IterativeLSI3",
    method="Seurat",
    name="Clusters6",
    resolution=0.8,
    force=TRUE
)



saveArchRProject(
    ArchRProj=proj2,
    outputDirectory="Gonad_ArchRProject-beforeUMAP",
    overwrite=TRUE,
    load=T,
    dropCells = FALSE,
    logFile=createLogFile("saveArchRProject"),
    threads=10
)

proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 20,
    minDist = 0.8,
    metric = "cosine",
    force=TRUE
)
proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "IterativeLSI2",
    name = "UMAP2",
    nNeighbors = 20,
    minDist = 0.8,
    metric = "cosine",
    force=TRUE
)
proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "IterativeLSI3",
    name = "UMAP3",
    nNeighbors = 20,
    minDist = 0.8,
    metric = "cosine",
    force=TRUE
)

proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 20,
    minDist = 0.8,
    metric = "cosine",
    force=TRUE
)
proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "Harmony2",
    name = "UMAPHarmony2",
    nNeighbors = 30,
    minDist = 0.4,
    metric = "cosine",
    force=TRUE
)
proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "Harmony3",
    name = "UMAPHarmony3",
    nNeighbors = 20,
    minDist = 0.8,
    metric = "cosine",
    force=TRUE
)

saveArchRProject(
    ArchRProj=proj2,
    outputDirectory="Gonad_ArchRProject-afterUMAP",
    overwrite=TRUE,
    load=T,
    dropCells = FALSE,
    logFile=createLogFile("saveArchRProject"),
    threads=12
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "gender", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "stage", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters3", embedding = "UMAP")
#p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,name="Plot-UMAP-PeakMatrix.0205.pdf",ArchRProj = proj, addDOC = FALSE, width = 50, height = 50)


proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI2",
    iterations = 4,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    seed = 1,
    force = TRUE
)
proj <- addHarmony(
    ArchRProj=proj,
    reducedDims="IterativeLSI2",
    name="Harmony2",
    groupBy="stage",
    force=TRUE
)
proj <- addClusters(
    input=proj,
    reducedDims="Harmony2",
    method="Seurat",
    name="Clusters4",
    resolution=0.8,
    force=T)
proj3 <- addUMAP(
    ArchRProj=proj3,
    reducedDims="IterativeLSI2",
    name="UMAP2final",
    nNeighbors=20,
    minDist=0.6,
    metric="cosine",
    force=TRUE)
proj3 <- addUMAP(
    ArchRProj=proj3,
    reducedDims="Harmony2",
    name="UMAPHarmonyfinal",
    nNeighbors=50,
    minDist=2.0,
    metric="cosine",
    force=TRUE)


p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony3")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony3")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony3")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters4", embedding = "UMAPHarmony3")
#p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,name="Plot-UMAPHarmony2-TileMatrix.0205.pdf",ArchRProj = proj, addDOC = FALSE, width = 50, height = 50)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmonyfinal")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmonyfinal")
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmonyfinal")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters5", embedding = "UMAPHarmonyfinal")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters7", embedding = "UMAPHarmonyfinal")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAPHarmony-PeakMatrix.0207.FINAL.pdf",ArchRProj = proj3, addDOC = FALSE, width = 50, height = 50)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAP3")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "gender", embedding = "UMAP3")
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "stage", embedding = "UMAP3")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters6", embedding = "UMAP3")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP3")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAP3-GeneScoreMatrix.0202.pdf",ArchRProj = proj3, addDOC = FALSE, width = 50, height = 50)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAPHarmony-TileMatrix.0202.pdf",ArchRProj = proj3, addDOC = FALSE, width = 50, height = 50)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony2")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony2")
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony2")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony2")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAPHarmony2")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAPHarmony2-PeakMatrix.0202.2.pdf",ArchRProj = proj3, addDOC = FALSE, width = 50, height = 50)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony3")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony3")
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony3")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters3", embedding = "UMAPHarmony3")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAPHarmony3")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAPHarmony3-GeneScoreMatrix.0202.pdf",ArchRProj = proj3, addDOC = FALSE, width = 50, height = 50)


## Identifying Marker Features
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup_Co",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerPeakList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR=TRUE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
  transpose = TRUE
)

#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot",annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peaks-Marker-Heatmap",width=10,height=8,ArchRProj=proj3, addDOC=FALSE)

proj3$Clusters7 <- proj3$predictedGroup_Co
proj3$Clusters7 <- gsub("Pre_Granulosa_Cell","01_Pre_Granulosa_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Sertoli_Cell","02_Sertoli_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Bipotential_State","03_Bipotential_State",proj3$Clusters7)
proj3$Clusters7 <- gsub("Fetal_Leydig_Cell","06_Fetal_Leydig_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Mesothelial_Cell","05_Mesothelial_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Stromal_Progenitor_Cell","04_Stromal_Progenitor_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Primordial_Germ_Cell","07_Primordial_Germ_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Meiotic_Germ_Cell","08_Meiotic_Germ_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Endothelial_Cell","09_Endothelial_Cell",proj3$Clusters7)
proj3$Clusters7 <- gsub("Erythroblast","10_Erythroblast",proj3$Clusters7)
proj3$Clusters7 <- gsub("Macrophage","11_Macrophage",proj3$Clusters7)

#Plot-Tracks-marker-peaks
p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Clusters7", 
    #geneSymbol = c("Kctd14","Lgals7","Fst","Amh","Clec18a","Col1a2","Upk3b","Cyp17a1","Dppa5a","Tex12","Cd34","Hbb-bs","C1qb"),
    #geneSymbol=c("Fst","Sox9","Myl7","Krt7"),
    geneSymbol=c("Polr3a","Cebpb","Smyd3"),
    features =  markerPeakList,
    upstream = 300000,
    downstream = 300000
)
plotPDF(p, name = "Plot-Tracks-With-Features-Mesothelial.pdf", width = 12, height = 10, ArchRProj = proj3, addDOC = FALSE)

#GRanges to dataframe
for (i in rownames(table(proj3$predictedGroup_Co))){
gr <- markerPeakList@listData[i]
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  strands=strand(gr),
  Log2FC=elementMetadata(gr)$Log2FC,
  FDR=elementMetadata(gr)$FDR,
  MeanDiff=elementMetadata(gr)$MeanDiff
)
write.table(df, file=paste0(i,".markerPeak.tsv"),quote=F,sep="\t",row.names=F, col.names=T)
}

#Trajectory
slabel <- rownames(table(proj3$Clusters7))
sclabel <- slabel[grep("01|02|03",slabel)]
idxSample <- BiocGenerics::which(proj3$Clusters7 %in% sclabel)
cellSample <- proj3$cellNames[idxSample]
proj_sc <- proj3[cellSample,]

proj_sc <- addIterativeLSI(
    ArchRProj = proj_sc,
    useMatrix = "PeakMatrix",
    name = "IterativeLSI_sc",
    iterations = 3,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    seed = 2
)
proj_sc <- addHarmony(
    ArchRProj=proj_sc,
    reducedDims="IterativeLSI_sc",
    corCutOff=0.7,
    name="Harmony_sc",
    groupBy="stage",
    force=TRUE
)
#Run UMAP
proj_sc <- addUMAP(
    ArchRProj = proj_sc,
    reducedDims = "IterativeLSI_sc",
    name = "UMAP_sc", 
    nNeighbors = 60,
    minDist = 2.1,
    metric = "cosine",
    force = TRUE,
    seed=20,
    threads=10
)
p1 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "sample", embedding = "UMAP_sc")
p2 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "Clusters7", embedding = "UMAP_sc")
p3 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "gender", embedding = "UMAP_sc")
p4 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "stage", embedding = "UMAP_sc")
plotPDF(p1,p2,p3,p4,p5,name="Plot-UMAP-SupportingCells.0207.6.pdf", ArchRProj = proj_sc, addDOC = FALSE, width = 20, height = 25)

proj_sc <- addUMAP(
    ArchRProj = proj_sc,
    reducedDims = "Harmony_sc",
    name = "UMAPHarmony_sc",
    nNeighbors = 50,
    minDist = 1.5,
    metric = "cosine",
    force=TRUE
)
p1 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "sample", embedding = "UMAPHarmony_sc")
p2 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "Clusters7", embedding = "UMAPHarmony_sc")
p3 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "gender", embedding = "UMAPHarmony_sc")
p4 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "stage", embedding = "UMAPHarmony_sc")
plotPDF(p1,p2,p3,p4,name="Plot-UMAP2Harmony-SupportingCells.0207.pdf", ArchRProj = proj_sc, addDOC = FALSE, width = 5, height = 5)

proj_sc <- addClusters(
    input = proj_sc,
    reducedDims = "IterativeLSI_sc",
    method = "Seurat",
    name = "Clusters_sc",
    resolution = 0.8,
    force = TRUE
)
p5 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "Clusters_sc", embedding = "UMAP_sc")

proj_sc$ClusterG <- paste0(proj_sc$Clusters_sc,proj_sc$gender)
table(proj_sc$Clusters_sc,proj_sc$gender)
#        XX   XY
#  C1    49  311
#  C10 1931   56
#  C11  760   30
#  C12    2   99
#  C13   12  179
#  C14  268  174
#  C2   592  249
#  C3   203  111
#  C4    66  983
#  C5     0 1547
#  C6   817  203
#  C7  2449   36
#  C8  1944   34
#  C9  2510    2

#New Groups
proj_sc$ClusterNEW <- paste0(proj_sc$Clusters7,proj_sc$gender)
GroupS <- c("03_Pre_Supporting_CellXY","02_Sertoli_CellXX","02_Sertoli_CellXY")
GroupG <- c("03_Pre_Supporting_CellXX","01_Pre_Granulosa_CellXY","01_Pre_Granulosa_CellXX")

sce <- SingleCellExperiment(
    assays = SimpleList(
      counts = as(matrix(rnorm(nCells(proj_sc) * 3), ncol = nCells(proj_sc), nrow = 3), "dgCMatrix")
    ),
    colData = getCellColData(proj_sc)
  )

  cds <- methods::new(
    "cell_data_set", 
    assays = SummarizedExperiment::Assays(list(counts = methods::as(assay(sce), "dgCMatrix"))), 
    colData = colData(sce), 
    int_elementMetadata = int_elementMetadata(sce), 
      int_colData = int_colData(sce), 
      int_metadata = int_metadata(sce), 
      metadata = metadata(sce), 
      NAMES = NULL, 
      elementMetadata = elementMetadata(sce)[, 0], 
      rowRanges = rowRanges(sce)
  )

cds@reduce_dim_aux@listData[["UMAP"]] <- getEmbedding(ArchRProj, embedding = embedding)

  if(!is.null(useGroups)){
    cds <- cds[, which(colData(cds)[, groupBy] %in% useGroups)]
  }



proj_sc2 <- getMonocleTrajectories(
  ArchRProj = proj_sc,
  name = "MonoTrajS",
  useGroups = GroupS,
  principalGroup = "03_Pre_Supporting_CellXY",
  groupBy = "ClusterNEW",
  embedding = "UMAPHarmony_sc",
  clusterParams = list(),
  graphParams = list(),
  seed = 1
)

> proj_sc$ClusterT <- proj_sc$ClusterG
> proj_sc$ClusterT <- gsub("C2XX","S1XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C3XX","S1XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C1XX","S2XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C12XX","S3XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C13XX","S4XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C2XY","S1XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C3XY","S1XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C1XY","S2XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C12XY","S3XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C13XY","S4XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C14XY","S5XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C14XX","S5XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C7XX","S6XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C7XY","S6XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C8XY","S6XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C8XX","S6XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C9XX","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C9XY","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C10XY","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C10XX","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C11XX","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C11XY","S7XX",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C4XY","S5XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C4XX","S5XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C5XX","S6XY",proj_sc$ClusterT)
> proj_sc$ClusterT <- gsub("C5XY","S6XY",proj_sc$ClusterT)
> table(proj_sc$ClusterT)

p6 <- plotEmbedding(ArchRProj = proj_sc, colorBy = "cellColData", name = "ClusterT", embedding = "UMAP_sc")


#Sertoli
traj1 <- c("S1XY","S2XY","S3XY","S4XY","S5XY","S6XY")

#Granulosa
traj2 <- c("S1XX","S2XX","S3XX","S4XX","S5XX","S6XX","S7XX")

proj_sc <- addTrajectory(
    ArchRProj = proj_sc,
    name = "Bipotential2SertoliU",
    groupBy = "ClusterT",
    trajectory = traj1,
    embedding = "UMAP_sc",
    force = TRUE
)
p1 <- plotTrajectory(proj_sc, trajectory = "Bipotential2SertoliU", colorBy = "cellColData", name = "Bipotential2SertoliU",embedding="UMAP_sc")

proj_sc <- addTrajectory(
    ArchRProj = proj_sc,
    name = "Bipotential2GranulosaU",
    groupBy = "ClusterT",
    trajectory = traj2,
    embedding = "UMAP_sc",
    force = TRUE
)
p2 <- plotTrajectory(proj_sc, trajectory = "Bipotential2GranulosaU", colorBy = "cellColData", name = "Bipotential2GranulosaU",embedding="UMAP_sc")

plotPDF(p1,p2, name = "Plot-Traj-UMAP.0207.pdf", ArchRProj = proj_sc, addDOC = FALSE, width = 10, height = 10)

proj_sc$SampleNum <- proj_sc$stage
proj_sc$SampleNum <- gsub("E11.5",1,proj_sc$SampleNum)
proj_sc$SampleNum <- gsub("E12.5",2,proj_sc$SampleNum)
proj_sc$SampleNum <- gsub("E13.5",3,proj_sc$SampleNum)
proj_sc$SampleNum <- as.numeric(proj_sc$SampleNum)

markerGenes <- c("Amh","Cst9","Sox9","Ptgds","Fst","Irx3","Gng13")
for (i in markerGenes){
p1 <- plotTrajectory(proj_sc, trajectory = "Bipotential2SertoliU",
        colorBy = "GeneScoreMatrix", name = i,
        continuousSet = "horizonExtra", embedding = "UMAP_sc")
p2 <- plotTrajectory(proj_sc, trajectory = "Bipotential2SertoliU",
        colorBy = "GeneIntegrationMatrix", name = i,
        continuousSet = "blueYellow", embedding = "UMAP_sc")
#ggAlignPlots(p1[[1]], p2[[1]], type = "h")
#ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1[[1]],p2[[1]],p1[[2]],p2[[2]], name = paste0("Trajectory-SertoliU-",i,".1202.pdf"), ArchRProj = proj_sc,
    addDOC = FALSE, width = 5, height = 5)
}

for (i in markerGenes){
p1 <- plotTrajectory(proj_sc, trajectory = "Bipotential2GranulosaU",
        colorBy = "GeneScoreMatrix", name = i,
        continuousSet = "horizonExtra", embedding = "UMAP_sc")
p2 <- plotTrajectory(proj_sc, trajectory = "Bipotential2GranulosaU",
        colorBy = "GeneIntegrationMatrix", name = i,
        continuousSet = "blueYellow", embedding = "UMAP_sc")
#ggAlignPlots(p1[[1]], p2[[1]], type = "h")
#ggAlignPlots(p1[[2]], p2[[2]], type = "h")
plotPDF(p1[[1]],p2[[1]],p1[[2]],p2[[2]], name = paste0("Trajectory-GranulosaU-",i,".0207.pdf"), ArchRProj = proj_sc,
    addDOC = FALSE, width = 5, height = 5)
}




p1 <- plotEmbedding(
    ArchRProj = proj_sc, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP_sc",
    imputeWeights = getImputeWeights(proj_sc)
)
plotPDF(plotList = p1,
    name = "Plot-UMAPHarmony-Marker-Genes-WO-Imputation-1.0207.pdf",
    ArchRProj = proj_sc,
    addDOC = FALSE, width = 5, height = 5)

a <- getMatrixFromProject(
  ArchRProj = proj_sc,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = 12,
  logFile = createLogFile("getMatrixFromProject")
)
