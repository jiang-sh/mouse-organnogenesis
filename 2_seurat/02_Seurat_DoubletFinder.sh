#!/bin/bash

Sample=HindBrain_E13.5
ratio=0.1

Organ=${Sample%_*}
Stage=${Sample#*_}
outputid=${Sample/./_}

rscript=/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/Scripts/dsub_02_Seurat_DoubletFinder/DoubletFinder_${outputid}.R

cat >${rscript} <<EOF
library('Seurat')
library('dplyr')
library('DoubletFinder')
library('KernSmooth')

scRNAin = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_Seurat_unannotated_rds/Seurat_scRNA_${outputid}_res_0.6.0323.rds"
scRNAout = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_Seurat_rmDoublet_rds/Seurat_scRNA_${outputid}.rmDoublet.rds"
sample <- readRDS(scRNAin)

# DoubletFinder
#pK Identification(no ground-truth)
sweep.res.list <- paramSweep_v3(sample, PCs=1:10, sct=FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
bcmvn <- find.pK(sweep.stats)
#Homotypic Doublet Proportion Estimate
annotations <- sample@meta.data\$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(${ratio}*length(annotations))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#Run DoubletFinder with varying calssification stringencies
sample <- doubletFinder_v3(sample, PCs=1:10, pN=0.25, pK=0.09, nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)

#***************************************************************
#*******Check whether pANN is same******************************
#***************************************************************

sample <- doubletFinder_v3(sample, PCs=1:10, pN=0.25, pK=0.09, nExp=nExp_poi.adj, reuse.pANN=colnames(sample@meta.data)[16], sct=FALSE)
sample@meta.data[,"DF_results"] <- sample@meta.data[17]
sample@meta.data\$DF_results[which(sample@meta.data\$DF_results == "Doublet" & sample@meta.data[18] == "Singlet")] <- "Doublet_low"
sample@meta.data\$DF_results[which(sample@meta.data\$DF_results == "Doublet")] <- "Doublet_high"

pdf("/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/03_Seurat_rmDoublet_rds/figure_pdf/DoubletFinder_${outputid}.pdf", width=10, height=5)
plot1 <- DimPlot(sample, reduction="umap", label=TRUE)
plot2 <- DimPlot(sample, group.by="DF_results", reduction="umap", plot.order=c("Doublet_high","Doublet_low","Singlet"), colors.use=c("black","gold","red"))
CombinePlots(plots=list(plot1,plot2), legend='bottom')
dev.off()
table(sample@meta.data\$DF_results)
aaa <- levels(sample\$seurat_clusters)
for (i in aaa){
a <- table(sample@meta.data\$DF_results[which(sample@meta.data\$seurat_clusters==i)])['Doublet_high']
b <- length(sample@meta.data\$DF_results[which(sample@meta.data\$seurat_clusters==i)])
print(a/b)
}

#*******************************************************************
#**************Check the cluster need to remove*********************
#*******************************************************************

#Remove Doublets cluster
sample <- subset(sample, subset = seurat_clusters!=11)

#Remove Gm42418_high cluster
sample <- subset(sample, subset = seurat_clusters!=7)
sample <- subset(sample, subset = seurat_clusters!=12)

saveRDS(sample, file=scRNAout)

EOF

chmod 755 ${rscript}

