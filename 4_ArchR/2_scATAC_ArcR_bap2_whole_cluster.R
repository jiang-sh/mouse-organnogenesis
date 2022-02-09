library(ArchR)
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/whole_Tss12000_v3")

addArchRThreads(threads = 8)
addArchRGenome("mm10")

fragments.files = c()

sample.names = c(
  "Lung_E12_5", "Liver_E13_5", "Pancreas_E13_5", "Spleen_E13_5", "ForeBrain_E13_5",
  "MidBrain_E13_5", "HindBrain_E13_5", "SpinalCord_E13_5", "Lung_E13_5", "Stomach_E13_5",
  "Intestine_E13_5", "Kidney_E13_5", "Gonad-male_E13_5", "Gonad-female_E13_5", "ForeLimb_E13_5",
  "Gonad-male_E12_5", "Gonad-female_E12_5", "Heart_E13_5", "Eye_E13_5", "GermLayer_E7_5",
  "Gonad-male_E11_5", "Gonad-female_E11_5"
);


for(index in 1:length(sample.names)){
  name_i = sample.names[index]
  fragments.files[index] = gsub("organ",name_i,"/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_barcode_multiplet_ArchR_bap2_v2/scATAC_organ_bap/final/scATAC_organ.fragments.tsv.gz")
}

## 1_create arrow file
ArrowFiles <- createArrowFiles(
  inputFiles = fragments.files,
  sampleNames = paste0("scATAC_", sample.names),
  filterTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = 0,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

## check arrow file
ArrowFiles

## 2_Doublet finding
#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#  LSIMethod = 1
#)

#print("finish doublet finding")


# 3.1_Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "ArchRProject",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

table(proj$Sample)

### use filter result from each sample ##
cellid <- read.table("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/08_scATAC_ArchR/sample_ArchR_bap2_fragment_v3/sample_barcode/wholemap_QC_rmD_barcode.txt", comment.char="")
proj <- proj[cellid$V1,]
table(proj$Sample)

## ask which data matrices are available within the ArchRProject
getAvailableMatrices(proj)

## add columns to cellColData
organ <- gsub("_E13_5|_E12_5|_E11_5|_E7_5|scATAC_","",proj$Sample)
proj$organ <- organ

id = strsplit(proj$cellNames, '[#]')
cellid = unlist(lapply(seq(id), function(x){id[[x]][2]}))
proj$cellid <- cellid

## add batch to whole map
batch = proj$Sample
batch = ifelse(batch == "scATAC_Lung_E12_5", "batch1", batch)
batch = ifelse(batch == "scATAC_Liver_E13_5", "batch2", batch)
batch = ifelse(batch == "scATAC_Pancreas_E13_5", "batch2", batch)
batch = ifelse(batch == "scATAC_Spleen_E13_5", "batch2", batch)
batch = ifelse(batch == "scATAC_ForeBrain_E13_5", "batch3", batch)
batch = ifelse(batch == "scATAC_MidBrain_E13_5", "batch3", batch)
batch = ifelse(batch == "scATAC_HindBrain_E13_5", "batch3", batch)
batch = ifelse(batch == "scATAC_SpinalCord_E13_5", "batch3", batch)
batch = ifelse(batch == "scATAC_Lung_E13_5", "batch4", batch)
batch = ifelse(batch == "scATAC_Stomach_E13_5", "batch4", batch)
batch = ifelse(batch == "scATAC_Intestine_E13_5", "batch4", batch)
batch = ifelse(batch == "scATAC_Kidney_E13_5", "batch4", batch)
batch = ifelse(batch == "scATAC_Gonad-female_E13_5", "batch5", batch)
batch = ifelse(batch == "scATAC_Gonad-male_E13_5", "batch5", batch)
batch = ifelse(batch == "scATAC_ForeLimb_E13_5", "batch6", batch)
batch = ifelse(batch == "scATAC_Gonad-male_E12_5", "batch7", batch)
batch = ifelse(batch == "scATAC_Gonad-female_E12_5", "batch7", batch)
batch = ifelse(batch == "scATAC_Heart_E13_5", "batch8", batch)
batch = ifelse(batch == "scATAC_Eye_E13_5", "batch8", batch)
batch = ifelse(batch == "scATAC_GermLayer_E7_5", "batch9", batch)
batch = ifelse(batch == "scATAC_Gonad-male_E11_5", "batch10", batch)
batch = ifelse(batch == "scATAC_Gonad-female_E11_5", "batch10", batch)
proj$batch <- batch

table(proj$organ)
table(proj$batch)

saveArchRProject(ArchRProj = proj, outputDirectory = "Save_scATAC_whole_ArchR_bap2_QC", load = FALSE)

# 3.2_filter doublets 
#proj2 <- filterDoublets(proj)
#saveArchRProject(ArchRProj = proj2, outputDirectory = "Save_Proj_rm_doublet", load = FALSE)

proj2 <- proj

## grep cell barcode for each sample after doublet filter
#ArchRSample <- paste0("scATAC_", sample.names)
#for (i in 1:length(ArchRSample)){
#  idxSample <- BiocGenerics::which(proj2$Sample %in% ArchRSample[i])
#  cellsSample <- proj2$cellid[idxSample]
#  write.table(cellsSample, file = paste0("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/08_scATAC_ArchR/scATAC_ArchR_bap2_fragment/cellsSampleBarcode/rm_doublet_",ArchRSample[i],".txt"), quote = F, col.names = F, row.names = F)
#}

## subset the project to keep all cells corresponding to a specific sample
#for (i in 1:length(ArchRSample)){
#  idxSample <- BiocGenerics::which(proj2$Sample %in% ArchRSample[i])
#  cellsSample <- proj2$cellNames[idxSample]
#  sampleproj <- proj2[cellsSample, ]
#  saveArchRProject(ArchRProj = proj2, outputDirectory = paste(ArchRSample[i],"_rm_doublet"), load = FALSE)
#}

# 3.3_Plotting Sample Statistics from an ArchRProject
## 3.3.1 Make a ridge plot for each sample for the TSS enrichment scores
p1 <- plotGroups(
    ArchRProj = proj2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

##  Make a violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
    ArchRProj = proj2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
## 3.3.2 Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
    ArchRProj = proj2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
## Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
    ArchRProj = proj2, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics-whole.pdf", ArchRProj = proj2, addDOC = FALSE, width = 4, height = 4)

## 3.3.3 Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles
p5 <- plotFragmentSizes(ArchRProj = proj2)

p6 <- plotTSSEnrichment(ArchRProj = proj2)
plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile-whole.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

save.image("scATAC_ArchR_bap2_fragment_before_dim.RData")

# 4_Dimensionality Reduction with ArchR  ## use projection to speed up 
## name parameter correspond to reducedDims object
proj2 <- addIterativeLSI(
    ArchRProj = proj2,
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
    sampleCellsPre = 30000,
    projectCellsPre = TRUE
)

## Batch Effect Correction wtih Harmony
## harmony batch and organ
proj3 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony2",
    groupBy = c("batch","organ")
)

# 5_Clustering with ArchR
proj3 <- addClusters(
    input = proj3,
    reducedDims = "Harmony2",
    method = "Seurat",
    name = "Clusters",
    resolution = 7
)

cM3 <- confusionMatrix(paste0(proj3$Clusters), paste0(proj3$Sample))

save.image("scATAC_ArchR_bap2_fragment_before_addUmap.RData")

library(pheatmap)
pdf("cluster_sample_distribution.pdf")
cM3 <- cM3 / Matrix::rowSums(cM3)
p3 <- pheatmap::pheatmap(
    mat = as.matrix(cM3),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)
p3

dev.off()

## 6_run UMAP
proj3 <- addUMAP(
    ArchRProj = proj3, 
    reducedDims = "Harmony2", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1,p2, name = "Plot-UMAPHarmony-Sample-Clusters-harmonyBatchOrgan.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj3, outputDirectory = "Save_scATAC_whole_ArchR_bap2_harmonyBatchOrgan", load = FALSE)
save.image("scATAC_ArchR_bap2_fragment_UMAP.RData")
