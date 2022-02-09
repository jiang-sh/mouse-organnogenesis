#Use R3.6

library("Seurat")

path = "/asnas/liujiang_group/yuchwei/Mouse_Organogenesis_v2/01_Seurat_wholeMap/subcluster/rds/Seurat_wholeMap_harmony_cluster"
#i=45
for (i in 0:45){
  pathid=paste0(path,i)
  print(pathid)
  sample <- readRDS(paste0(pathid,".rds"))
  cat(paste0(paste0("Cluster_",i),"\n"),file="DoubletRatio_subclusters.results",append=TRUE)
  print(pathid)
  table(sample@meta.data$DF_results)
  aaa <- levels(sample$seurat_clusters)
  for (j in aaa){
    a <- table(sample@meta.data$DF_results[which(sample@meta.data$seurat_clusters==j)])['Doublet_high']
    b <- length(sample@meta.data$DF_results[which(sample@meta.data$seurat_clusters==j)])
    cat(paste0(a/b,"\n"),file="DoubletRatio_subclusters.results",append=TRUE)
  }
}
