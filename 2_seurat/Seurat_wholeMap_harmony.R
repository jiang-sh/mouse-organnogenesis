library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(biomaRt)

pdf.file = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap_pdf/Seurat_wholeMap_harmony_cell"
rds = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/Seurat_wholeMap_harmony_cell"

sample.names = c(
  "Lung_E12_5", "Liver_E13_5", "Pancreas_E13_5", "Spleen_E13_5", "ForeBrain_E13_5", 
  "MidBrain_E13_5", "HindBrain_E13_5", "SpinalCord_E13_5", "Lung_E13_5", "Stomach_E13_5", 
  "Intestine_E13_5", "Kidney_E13_5", "Gonad-male_E13_5", "Gonad-female_E13_5", "ForeLimb_E13_5",
  "Gonad-male_E12_5", "Gonad-female_E12_5", "Heart_E13_5", "Eye_E13_5", "GermLayer_E7_5",  
  "Gonad-male_E11_5", "Gonad-female_E11_5"
);
germlayer = c("Endoderm", "Endoderm", "Endoderm", "Mesoderm", "Ectoderm",
              "Ectoderm", "Ectoderm", "Ectoderm", "Endoderm", "Endoderm",
              "Endoderm", "Mesoderm", "Mesoderm", "Mesoderm", "Mesoderm",
              "Mesoderm", "Mesoderm", "Mesoderm", "Ectoderm", "E7.5",
              "Mesoderm", "Mesoderm")
sample.batch = c(
  "batch1", "batch2", "batch2", "batch2", "batch3", 
  "batch3", "batch3", "batch3", "batch4", "batch4", 
  "batch4", "batch4", "batch5", "batch5", "batch6",
  "batch7",  "batch7", "batch8", "batch8", "batch9", 
  "batch10", "batch10"
);
cell.id = c(
  "E12-LU","E13-LI","E13-PA","E13-SP","E13-FB",
  "E13-MB","E13-HB","E13-SC","E13-LU","E13-ST",
  "E13-IN","E13-KI","E13-GM","E13-GF","E13-FL",
  "E12-GM","E12-GF","E13-HE","E13-EY","E7-GL",
  "E11-GM","G11-GF"
)

#### read filtered doublet & unannotated rna rds, add batch.id and sample name to seurat object####
print("read rds file")
rds.files = c()
rds.ls = list()
for(i in 1:length(sample.names)){
  name_i = sample.names[i]
  rds.files[i] = gsub("organ",name_i,"/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_Seurat_rmDoublet_rds/Seurat_scRNA_organ.rmDoublet.rds")
  rds.ls[[cell.id[i]]] <- readRDS(file = rds.files[i])
  rds.ls[[cell.id[i]]]@meta.data$batch = sample.batch[i]
  rds.ls[[cell.id[i]]]@meta.data$sample = cell.id[i]
  rds.ls[[cell.id[i]]]@meta.data$germlayer = germlayer[i]
}
rds.ls

# merge all seurat object
print("merge seurat rds files")
all <- merge(x = rds.ls$`E12-LU`, 
             y = c(rds.ls$`E13-LI`, rds.ls$`E13-PA`, rds.ls$`E13-SP`, rds.ls$`E13-FB`, 
                   rds.ls$`E13-MB`, rds.ls$`E13-HB`, rds.ls$`E13-SC`, rds.ls$`E13-LU`, rds.ls$`E13-ST`,
                   rds.ls$`E13-IN`, rds.ls$`E13-KI`, rds.ls$`E13-GM`, rds.ls$`E13-GF`, rds.ls$`E13-FL`,
                   rds.ls$`E12-GM`, rds.ls$`E12-GF`, rds.ls$`E13-HE`, rds.ls$`E13-EY`, rds.ls$`E7-GL`,
                   rds.ls$`E11-GM`, rds.ls$`G11-GF`), 
             add.cell.ids = cell.id, 
             merge.data = T, 
             project = "Mouse atlas")

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- read.table("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/cellcycle_s.txt")
g2m.genes <- read.table("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/00_reference/cellcycle_g2m.txt")

# 1.Initialize the Seurat object with the raw (non-normalized data) 
all

########## Standard pre-processing workflow ###########
# 1_QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
print("calculate mt percentage")
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")       #all@meta.data$percent.mt

# Visualize QC metrics as a violin plot
pdf(paste0(pdf.file, "_01_QC.pdf"))
options(repr.plot.height = 2, repr.plot.width = 18)
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")

options(repr.plot.height = 4, repr.plot.width = 6)
#UMI, gene, per.mt distribution
hist(log10(all@meta.data$nCount_RNA), xlab = "UMIS per cell (log10)", ylab = "Frequency", main = "UMI distribution", breaks = 80, col = "mediumseagreen")
hist(log10(all@meta.data$nFeature_RNA), xlab = "Genes per cell (log10)", ylab = "Frequency", main = "Gene distribution", breaks = 80, col = "mediumseagreen")
hist(all@meta.data$percent.mt, xlim = c(0,20),xlab = "Per.mt per cell", ylab = "Frequency", main = "Per.mt distribution", breaks = 180, col = "mediumseagreen")

options(repr.plot.height = 2, repr.plot.width = 12)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#all <- subset(all, subset = nFeature_RNA>500 & percent.mt<10)
print("cell number after quality control")
all
table(all@meta.data[["sample"]])

########## 2_Normalizing the data ########################
# “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
print("Normalizing the data")
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)

########## 3_Identification of highly variable features (feature selection)#####################
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(all), 10)

# plot variable features with and without labels
pdf(paste0(pdf.file, "_02_VariableFeature.pdf"))
options(repr.plot.height = 2, repr.plot.width = 12)
plot1 <- VariableFeaturePlot(all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

######### 4_Scaling the data ############
print("scale data")
##Assign Cell-Cycle Scores
all <- CellCycleScoring(all, s.features = s.genes$V1, g2m.features = g2m.genes$V1, set.ident = TRUE)
head(all[[]])
#cell-cycle scoring and regression
all.genes <- rownames(all)
all$CC.Difference <- all$S.Score - all$G2M.Score
all <- ScaleData(all, features = all.genes, vars.to.regress = c("CC.Difference"))

########## 5_Perform linear dimensional reduction ################################
print("perform PCA")
all <- RunPCA(all, features = VariableFeatures(object = all))

# Examine and visualize PCA results a few different ways
pdf(paste0(pdf.file, "_03_VisualizePCA.pdf"))
print(all[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all, dims = 1:2, reduction = "pca")
DimPlot(all, reduction = "pca", group.by = "sample")
DimPlot(all, reduction = "pca", group.by = "germlayer")
DimPlot(all, reduction = "pca", group.by = "batch")
DimHeatmap(all, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

################ 6_Determine statistically significant principal components ##########
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#2
#all <- JackStraw(all, num.replicate = 100)
#all <- ScoreJackStraw(all, dims = 1:20)
#JackStrawPlot(all, dims = 1:15)
#3
pdf(paste0(pdf.file, "_04_ElbowPlot.pdf"))
ElbowPlot(all, ndims = 50)
dev.off()

saveRDS(all, file = paste0(rds, "_PCA.rds"))


############### 7_Run harmony ###################################
options(repr.plot.height = 2.5, repr.plot.width = 6)
pdf(paste0(pdf.file, "_05_harmony_reduction.pdf"))
all <- all %>% 
    RunHarmony("batch", plot_convergence = TRUE)

### harmony embeddings ###
harmony_embeddings <- Embeddings(all, 'harmony')
harmony_embeddings[1:5, 1:5]

DimPlot(object = all, reduction = "harmony", pt.size = .1, group.by = "batch", do.return = TRUE)
DimPlot(all, reduction = "harmony", group.by = "sample", pt.size = .1)
DimPlot(all, reduction = "harmony", group.by = "germlayer", pt.size = .1)
dev.off()

########## visualization ########################### 
all6 <- all %>% 
    RunUMAP(reduction = "harmony", dims = 1:50) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
    FindClusters(resolution = 0.6) %>% 
    identity()

pdf(paste0(pdf.file, "_05_harmony_cluster.pdf"))
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(all6, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(all6, reduction = "umap", label = T, group.by = "sample")
DimPlot(all6, reduction = "umap", label = F, group.by = "sample")
DimPlot(all6, reduction = "umap", label = T, group.by = "germlayer")
DimPlot(all6, reduction = "umap", label = T, group.by = "Phase")
dev.off()

saveRDS(all, file = paste0(rds, "_runharmony.rds"))
saveRDS(all6, file = paste0(rds, "_runharmony_r0.6.rds"))

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(all6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers, file = paste0(rds, "_clusters.markers.xls"), sep = "\t", quote = F, col.names = T, row.names = T)

print("End!")
