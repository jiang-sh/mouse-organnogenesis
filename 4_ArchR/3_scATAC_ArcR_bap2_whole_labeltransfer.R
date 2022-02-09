library(ArchR)
library(Seurat)
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/whole_Tss12000_v3")

addArchRThreads(threads = 8)
addArchRGenome("mm10")

sample.names = c(
  "Lung_E12_5", "Liver_E13_5", "Pancreas_E13_5", "Spleen_E13_5", "ForeBrain_E13_5",
  "MidBrain_E13_5", "HindBrain_E13_5", "SpinalCord_E13_5", "Lung_E13_5", "Stomach_E13_5",
  "Intestine_E13_5", "Kidney_E13_5", "Gonad-male_E13_5", "Gonad-female_E13_5", "ForeLimb_E13_5",
  "Gonad-male_E12_5", "Gonad-female_E12_5", "Heart_E13_5", "Eye_E13_5", "GermLayer_E7_5",
  "Gonad-male_E11_5", "Gonad-female_E11_5"
);

load("scATAC_ArchR_bap2_fragment_UMAP.RData")

## 8_Defining Cluster Identity with scRNA-seq
## Cross-platform linkage of scATAC-seq cells with scRNA-seq cells
#seRNA <- readRDS("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/wholeMap/Seurat_wholeMap_harmony_cell_runharmony_r0.6_subCluster.labeled.v5.rds")
seRNA <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/sample_ArchR_bap2_fragment_v3/groupList/Seurat_wholeMap_harmony_cell_runharmony_r0.6_subCluster.labeled.v6.rds")
DefaultAssay(seRNA) <- "RNA"
seRNA

groupList <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/sample_ArchR_bap2_fragment_v3/groupList/grouplist.new.rna.barcode")

organ <- gsub("_E13_5|_E12_5|_E11_5|_E7_5|scATAC_","",proj3$Sample)
proj3$organ <- organ

#### 8. defining cluster identity with scRNA-seq
## 8.2 Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj3,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony2",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "Label43",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
proj3$Clusters2 <- proj3$predictedGroup
proj3 <- addImputeWeights(proj3)

saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj3_Co_addToarrow", load = FALSE)
save.image("scATAC_ArchR_bap2_fragment_UMAP_labeltransfer.RData")

pal <- paletteDiscrete(values = seRNA$Label43)

print("plot of 2nd LT")
p2 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup",
    pal = pal,
    embedding = "UMAPHarmony"
)
plotPDF(p2, name = "Plot-UMAP-Integration_Co.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

p4 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedScore",
    #pal = pal,
    embedding = "UMAPHarmony"
)
plotPDF(p3,p4, name = "8.2-Plot-UMAP-Integration_Score.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

## save second round predicted celltype plot
p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAPHarmony", baseSize = 1)
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "predictedScore", embedding = "UMAPHarmony", baseSize = 1)
plotPDF(p1,p2, name = "8.2-Plot-UMAP-2nd-predicted-Clusters.pdf", ArchRProj = proj3, addDOC = FALSE, width = 10, height = 10)

print("histogram")
## plot prediction score  histogram
pdf("8.2-whole_harmony2_2nd_r7_Clusters_predict_score.pdf")
hist(
    proj3$predictedScore,
    xlab="prediction score",
    col="lightblue",
    xlim=c(0, 1),
    main="scATAC prediction score"
  );
abline(v=0.4, col="red", lwd=2, lty=2);
table(proj3$predictedScore > 0.4);
dev.off()

### re-add group information
id = strsplit(proj3$predictedGroup, '[=]')
group = unlist(lapply(seq(id), function(x){id[[x]][1]}))
proj3$group <- group

print("8.2-violin plot")
## PLOT violin plot
Phen=proj3$group
score=proj3$predictedScore
df=data.frame(Phen,score)
agg=aggregate(df,by=list(df$Phen),FUN=mean)
agg=agg[order(-agg$score),]
level=agg$Group.1
colors=agg$score
print(colors)
df$Phen=factor(df$Phen,levels=level)

p<-ggplot(data=df, aes(x=Phen,y=score,group=Phen))+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 60, hjust = 1, size = 7))+
labs(x="Cell Type",y="predicted score",title="Gonad-female_E11_5")+guides(fill=F)
#pb=p+geom_boxplot(aes(fill=Phen))
pv=p+geom_violin(scale='width',aes(fill=Phen))

ggsave(pv,file= "8.2-constrain_Cluster_predictedscore.pdf")

print("END-label transfer")


### 9. pseudo-bulk replicates in ArchR  #### proj3 to proj4
print("9.addGroupCoverages")
proj4 <- addGroupCoverages(ArchRProj = proj3, groupBy = "Clusters2")

### 10. calling peaks with ArchR
## 10.2 caliing peaks macs2 #### previous defined group and repilcate
#pathToMacs2 <- findMacs2()
print("10.addReproduciblePeakSet")
proj4 <- addReproduciblePeakSet(
    ArchRProj = proj4,
    groupBy = "Clusters2",
    pathToMacs2 = "/asnas/liujiang_group/yuhao/Software/anaconda3/bin/macs2"
)
getPeakSet(proj4)

print("10.Save peakcslling")
saveArchRProject(ArchRProj = proj4, outputDirectory = "Save-Proj4_peakcalling", load = FALSE)

## 10.4 add peak matrix  #### proj4 to proj5
print("10.4-addPeakMatrix")
proj5 <- addPeakMatrix(proj4)
getAvailableMatrices(proj5)

### 11.1  Identifying Marker Peaks with ArchR
print("11.1-getMarkerFeatures")
table(proj5$Clusters2)
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj5,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
##11.2 Plotting Marker Peaks in ArchR
#Marker Peak Heatmaps
print("11.2-markerHeatmap")
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "11.2-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)

save.image("11_peakcalling.RData")


