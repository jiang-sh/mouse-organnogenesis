library(Matrix)
library(ComplexHeatmap)
library(ggplot2)
library(Seurat)
library(ArchR)

## 1 ##
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/00_yuhao_peak2gene_link")
CT2Peak_01.m = readRDS("celltype_peak_from_callpeakfile.m.binary.RDS")
CT2Peak_01.m.simple = as.matrix(CT2Peak_01.m)

# cluster_id_old_new_df <- read.table(file = "cluster_id_old_new_chw_df.txt")
# cluster_order <- cluster_id_old_new_df[order(cluster_id_old_new_df$order),"cluster"]
# CT2Peak_01.m.simple <- CT2Peak_01.m.simple[, as.character(cluster_order)]
peak_mt_order_final <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/Script_new/0_whole/NMF_clustering/peak_cluster_mt_order_final.rds")
cluster <- colnames(peak_mt_order_final);rm(peak_mt_order_final)
CT2Peak_01.m.simple <- CT2Peak_01.m.simple[, cluster]


## 2 ##
CT2Peak_01.m.simple.norm2 = sweep(CT2Peak_01.m.simple, 1, 
                                  apply(CT2Peak_01.m.simple, 1, function(i){sqrt(i %*% i)}),`/`)
km2 = kmeans(CT2Peak_01.m.simple.norm2,centers = 20,iter.max = 20,nstart = 2,algorithm = "Hartigan-Wong")
km.class2 = fitted(km2,method=c("classes"))
CT2Peak_01.m.simple.reorder2 = CT2Peak_01.m.simple[sort.int(km.class2,index.return = T)$ix,]

pdf("callpeakfile_binary_heatmap_mo_mannul.pdf",width = 50,height = 200)
Heatmap(CT2Peak_01.m.simple.reorder2,cluster_columns = F,cluster_rows = F,show_row_names = F,col = c("white","black"))
dev.off()
save(list = c("CT2Peak_01.m.simple.norm2","km2","km.class2","CT2Peak_01.m.simple.reorder2"),file = "km20_Total_mo.RData")
lapply(1:20, function(i){
  km2 = kmeans(CT2Peak_01.m.simple.norm2,centers = 20,iter.max = 20,nstart = 2,algorithm = "Hartigan-Wong")
  km.class2 = fitted(km2,method=c("classes"))
  CT2Peak_01.m.simple.reorder2 = CT2Peak_01.m.simple[sort.int(km.class2,index.return = T)$ix,]
  pdf(paste0("Out/callpeakfile_binary_heatmap_mo_rep",i,".pdf"),width = 50,height = 200)
  ht = Heatmap(CT2Peak_01.m.simple.reorder2,cluster_columns = F,cluster_rows = F,show_row_names = F,col = c("white","black"))
  print(ht)
  dev.off()
  save(list = c("km2","km.class2","CT2Peak_01.m.simple.reorder2"),file = paste0("Out/km20_Total_mo_rep",i,".RData"))
})

## 3 ##
km.class2 = factor(km.class2, levels = c(9,13,20,4,10,6,3,2,
                                         17,11,14,5,12,
                                         19,16,15,7,8,18,1))
CT2Peak_01.m.simple.reorder2 = CT2Peak_01.m.simple[sort.int(km.class2,index.return = T)$ix,]
row.average2 = apply(t(t(CT2Peak_01.m.simple.reorder2[1:158410,])*c(1:271)),1,function(vi){
  return(sum(vi)/sum(vi>0))
})
CT2Peak_01.m.simple.reorder2 = CT2Peak_01.m.simple.reorder2[c(sort.int(row.average2,index.return = T,decreasing = T)$ix,158411:739952),]


## 4 ##
class.ls = sort(km.class2)
levels(class.ls) = 1:20
km.class.reorder = data.frame("peak"=row.names(CT2Peak_01.m.simple.reorder2),
                              "kmeans_idx"=class.ls)

C19_average_19 = Reduce(cbind,lapply(2:20, function(i){
  mati = CT2Peak_01.m.simple.reorder2[which(km.class.reorder$kmeans_idx == i),]
  return(colMeans(mati))
}))
colnames(C19_average_19) = 2:20
pdf("kms2celltype_heatmap.pdf",width = 24,height = 40)
Heatmap(C19_average_19[271:1,],cluster_columns = F,cluster_rows = F,column_names_max_height = unit(12, "cm"),
        show_row_names = T,col = c("white","black"))
dev.off()
save.image("1021.RData")

# 1 # split cell types
split_cell_types = c(0,18,5,44,29,14,8,3,8,5,10,13,11,6,21,10,6,17,4,17,6,5,12)
split_cell_types2 = unlist(lapply(1:23, function(i){sum(split_cell_types[1:i])+1}))

# 2 # split peaks according to splited cell types
get_part <- function(i){
  which(order(c(split_cell_types2,i))==24) - 1
}
row.average2 = apply(t(t(CT2Peak_01.m.simple.reorder2[1:158410,])*c(271:1)),1,function(vi){
  return(sum(vi)/sum(vi>0))
})
row.average2.order = unlist(lapply(row.average2, get_part))

# 3 # split kms1
split_peaks = c(0,c(5401,263,18633,16495,5833,5626,3112,8332,10737,14182,16698,5682,5175,12894,3998,4538,8318,3191,4431,3552,725,594))
split_peaks2 = unlist(lapply(1:23, function(i){sum(split_peaks[1:i])+1}))
# [1]      1   5402   5665  24298  40793  46626  52252  55364  63696  74433
# [11]  88615 105313 110995 116170 129064 133062 137600 145918 149109 153540
# [21] 157092 157817 158411

split_master_peaks = c(0,158410,11547,72250,52157,21332,20913,15891,16091,21914,31945,22366,
                       17005,27810,51488,29601,8096,10036,36750,35082,79268)
split_master_peaks2 = unlist(lapply(1:21, function(i){sum(split_master_peaks[1:i])+1}))
# [1]      1 158411 169958 242208 294365 315697 336610 352501 368592 390506
# [11] 422451 444817 461822 489632 541120 570721 578817 588853 625603 660685
# [21] 739953

# 4 # reorder + heatmap
sample_all = function(i){sample(i,length(i))}
CT2Peak_01.m.simple.reorder3 = CT2Peak_01.m.simple.reorder2[c(158411:169957,sample_all(5402:5664),# 2
                                                              169958:242207,sample_all(1:5401),sample_all(5665:24297),# 3
                                                              242208:294364,sample_all(24298:40792),# 4
                                                              294365:315696,# 5
                                                              315697:336609,# 6
                                                              336610:352500,# 7
                                                              352501:368591,# 8
                                                              368592:390505,sample_all(40793:46625),# 9
                                                              390506:422450,sample_all(46626:52251),sample_all(145918:149108),# 10
                                                              sample_all(52252:55363),422451:444816,sample_all(55364:63695),sample_all(63696:74432),# 11
                                                              444817:461821,sample_all(74433:88614),# 12
                                                              461822:489631,sample_all(88615:105312),# 13
                                                              489632:541119,sample_all(105313:110994),# 14
                                                              sample_all(110995:116169),
                                                              541120:570720,sample_all(116170:129063),# 15
                                                              sample_all(129064:133061),
                                                              570721:578816,sample_all(133062:137599),# 16
                                                              sample_all(137600:145917),sample_all(149109:153539),
                                                              578817:588852,sample_all(153540:157091),# 17
                                                              588853:625602,sample_all(157092:157816),# 18
                                                              625603:660684,sample_all(157817:158410),# 19
                                                              660685:739952),]
CT2Peak_01.m.simple.reorder4 = CT2Peak_01.m.simple.reorder2[c(158411:169957,5402:5664,# 2
                                                              169958:242207,1:5401,5665:24297,# 3
                                                              242208:294364,24298:40792,# 4
                                                              294365:315696,# 5
                                                              315697:336609,# 6
                                                              336610:352500,# 7
                                                              352501:368591,# 8
                                                              368592:390505,40793:46625,# 9
                                                              390506:422450,46626:52251,145918:149108,# 10
                                                              52252:55363,
                                                              422451:444816,55364:63695,63696:74432,# 11
                                                              444817:461821,74433:88614,# 12
                                                              461822:489631,88615:105312,# 13
                                                              489632:541119,105313:110994,# 14
                                                              110995:116169,
                                                              541120:570720,116170:129063,# 15
                                                              129064:133061,
                                                              570721:578816,133062:137599,# 16
                                                              137600:145917,149109:153539,
                                                              578817:588852,153540:157091,# 17
                                                              588853:625602,157092:157816,# 18
                                                              625603:660684,157817:158410,# 19
                                                              660685:739952),]

pdf("callpeakfile_binary_heatmap_rm1.pdf",width = 50,height = 200)
Heatmap(CT2Peak_01.m.simple.reorder3,cluster_columns = F,cluster_rows = F,show_row_names = F,col = c("white","black"))
dev.off()

write.csv(data.frame("peak"=row.names(CT2Peak_01.m.simple.reorder3),
                     "kmeans_idx"=rep(1:23,c(11810,96284,68652,21332,20913,
                                             15891,16091,27747,40762,3112,41435,
                                             31187,44508,57170,5175,42495,
                                             3998,12634,12749,13588,37475,35676,79268))),
          file = "Data/peaks_order&kmeans.csv",row.names = F)
saveRDS(CT2Peak_01.m.simple.reorder3,file = "Data/CT2Peak_01.m.RDS")
CT2Peak_01.m.simple.reorder3 = readRDS(file = "Data/CT2Peak_01.m.RDS")

## 4 ##
# p2g.df_uniq_gene = read.table("peak_gene_link_selected_45213.txt",header = T,stringsAsFactors = F)
p2g.df_uniq_gene = read.table("k100_knn500_default/peak_gene_link_selected.txt",header = T,stringsAsFactors = F) # 48351
p2g.df_uniq_gene$peak = gsub(":|-","_",p2g.df_uniq_gene$peak)
peak_order = read.csv(file = "Data/peaks_order&kmeans.csv",stringsAsFactors = F)


p2g.df_uniq_gene$peak_km = peak_order$kmeans_idx[match(p2g.df_uniq_gene$peak,peak_order$peak)]
gene2km = tapply(p2g.df_uniq_gene$peak,list(p2g.df_uniq_gene$gene,p2g.df_uniq_gene$peak_km),length,default = 0)

gene2km.proportion = gene2km/apply(gene2km, 1, sum)
gene2km.proportion.df = as.data.frame(as.table(gene2km.proportion))

pdf("kms_proportion_of_genes.pdf",width = 15,height = 8)
ggplot(gene2km.proportion.df,aes(Var1,weight=Freq,fill=Var2))+geom_bar(position="stack")+labs(title="",y="Proportion")
dev.off()

gene_order = data.frame("gene"=row.names(gene2km),"Max_peak"=NA,"Max_count"=NA,"Top1_peak"=1,"Top1_km"=0,stringsAsFactors = F)
gene_order$Max_count = apply(gene2km,1,function(vi){
  v2 = sort(vi,decreasing = T)
  if(v2[2]==v2[1]){
    return(NA)
  }else{
    return(order(vi,decreasing = T)[1])
  }
})
gene_order$Max_peak[!is.na(gene_order$Max_count)] = unlist(lapply(which(!is.na(gene_order$Max_count)), function(i){
  p2g_i = p2g.df_uniq_gene[(p2g.df_uniq_gene$gene == gene_order$gene[i]) & (p2g.df_uniq_gene$peak_km == gene_order$Max_count[i]),]
  return(as.vector(p2g_i[which.max(p2g_i$Correlation),"peak"]))
}))

gene_order[,c("Top1_peak","Top1_km")] = Reduce(rbind,lapply(gene_order$gene, function(gi){
  p2g_i = p2g.df_uniq_gene[which(p2g.df_uniq_gene$gene == gi),]
  return(p2g_i[which.max(p2g_i$Correlation),c("peak","peak_km")])
}))
gene_order$Total_peak = gene_order$Max_peak
gene_order$Total_peak[is.na(gene_order$Max_peak)] = as.character(gene_order$Top1_peak[is.na(gene_order$Max_peak)])
gene_order$Total_km = gene_order$Max_count
gene_order$Total_km[is.na(gene_order$Total_km)] = gene_order$Top1_km[is.na(gene_order$Total_km)]
gene_order$Top1_peak = as.vector(gene_order$Top1_peak)
write.csv(gene_order,file = "Data/gene_order.csv")

gene_order$Total_peak = factor(gene_order$Total_peak,levels = as.vector(peak_order$peak))
table(gene_order$Total_peak %in% as.vector(peak_order$peak))
save(list = c("p2g.df_uniq_gene","peak_order","gene2km.proportion.df","gene_order"),
     file = "gene_order_originp2g.RData")

## 5 ##
cluster.average = readRDS("subcluster_gene_mt.rds")
table(gene_order$gene %in% colnames(cluster.average))
table(colnames(CT2Peak_01.m.simple.reorder3) %in% row.names(cluster.average))

rownames(cluster.average) <- gsub("[=]|[-]|[:]|[(]|[)]|[/]|[+]|[-]",".", rownames(cluster.average))
gene_order2 = gene_order
gene_order2$Total_peak = factor(gene_order2$Total_peak,levels = row.names(CT2Peak_01.m.simple.reorder3))
CT2gene.m2 = t(cluster.average)
CT2gene.m2 = CT2gene.m2[gene_order2$gene[order(gene_order2$Total_peak)],]
CT2gene.m2 = CT2gene.m2[,colnames((CT2Peak_01.m.simple.reorder3))]

CT2gene.m.scale2 = scale(as.matrix(t(CT2gene.m2)))
CT2gene.m.scale2[CT2gene.m.scale2 > 4] = 4
CT2gene.m.scale2[CT2gene.m.scale2 < -4] = -4
pdf("gene2celltype_heatmap_range8_zscore_cut4_originp2g.pdf",width = 200,height = 50)
Heatmap(CT2gene.m.scale2[271:1,],cluster_columns = F,cluster_rows = F,show_column_names = F,col = paletteContinuous("blueYellow"))
#Heatmap(CT2gene.m.scale2[271:1,],cluster_columns = F,cluster_rows = F,show_column_names = F,col = c("yellow","blue")) # forestgreen deepskyblue
dev.off()

## 6 ## ATAC
p2g.df_uniq_gene$gene = factor(p2g.df_uniq_gene$gene,levels = gene_order2$gene[order(gene_order2$Total_peak)])
CT2Peak_01.m.simple.reorder3.pairs = CT2Peak_01.m.simple.reorder3[p2g.df_uniq_gene$peak[order(p2g.df_uniq_gene$gene)],]
pdf("p2g_celltypepeak_heatmap_ATAC.pdf",width = 50,height = 200)
Heatmap(CT2Peak_01.m.simple.reorder3.pairs,cluster_columns = F,cluster_rows = F,show_row_names = F,col = c("white","black"))
dev.off()
saveRDS(CT2gene.m2,file = "Data/CT2gene.m.RData")
saveRDS(CT2Peak_01.m.simple.reorder3.pairs,file = "CT2Peak_01_pair.m.RDS")


## 7 ## RNA set max count
CT2gene.m2 = readRDS(file = "Data/CT2gene.m.RData")

pdf("gene2celltype_heatmap_range8_zscore_cut4_maxcount.pdf",width = 200,height = 50)
for (i in 1:5) {
  cluster.average.cut = readRDS(paste0("gene_expression_normalize/Out/",i,"00_average_expr.RDS"))
  cluster.average.cut = cluster.average.cut@assays$RNA@data
  colnames(cluster.average.cut) <- gsub("[=]|[-]|[:]|[(]|[)]|[/]|[+]|[-]",".", colnames(cluster.average.cut))
  cluster.average.cut = cluster.average.cut[row.names(CT2gene.m2),colnames(CT2gene.m2)]
  
  cluster.average.cut = scale(as.matrix(t(cluster.average.cut)))
  cluster.average.cut[cluster.average.cut > 4] = 4
  cluster.average.cut[cluster.average.cut < -4] = -4
  ht = Heatmap(cluster.average.cut[271:1,],cluster_columns = F,cluster_rows = F,show_column_names = F,
               col = paletteContinuous("blueYellow"),name = paste0("max count ",i,"00"))
  print(ht)
}
dev.off()


## 8 ##
library(stringr)
peak_order = read.csv("/home/yuhao/Desktop/Mouse/peaks_order&kmeans.csv",stringsAsFactors = F)
peak.ls = peak_order[which(peak_order$kmeans_idx==11),"peak"] # 9 11 18
peak.df = data.frame(str_split_fixed(peak.ls,'_',3),stringsAsFactors = F)
colnames(peak.df) = c("chr","start","end")
peak.df$name = peak.ls

peak.df$V5 = 1:length(peak.ls)
peak.df$V6 = "."
write.table(peak.df, file="/home/yuhao/Desktop/Mouse/kms11_col6.bed",
            row.names = F,col.names = F,sep = "\t",quote = F)
# findMotifsGenome.pl kms9_col6.bed mm10 km9 -size 500

lapply(1:23, function(i){
  peak.ls = peak_order[which(peak_order$kmeans_idx==i),"peak"] # 9 11 18
  peak.df = data.frame(str_split_fixed(peak.ls,'_',3),stringsAsFactors = F)
  colnames(peak.df) = c("chr","start","end")
  peak.df$name = peak.ls
  write.table(peak.df, file = paste0("Data/kms_peak/kms",i,"_col4.bed"),
              row.names = F,col.names = F,sep = "\t",quote = F)
})