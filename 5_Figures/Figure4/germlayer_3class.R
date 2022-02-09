############### 1. use binary matrix plot heatmap ##############
ct_mt <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/00_yuhao_peak2gene_link/celltype_peak_from_callpeakfile.m.binary.RDS")

germ_cluster_name <-read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/16_germlayer_organ/E7.5_germlayer_celltype.txt")
germ_cluster_name_new <- gsub("[=]|[-]|[:]|[(]|[)]|[/]|[+]|[-]",".", as.character(germ_cluster_name$V1))

three_germ_peak_uniq <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/germ_layer_specific_peak/three_germ_layer_specific_peak.bed")
three_germ_peak_uniq$V5 <- paste0(three_germ_peak_uniq$V1,"_",three_germ_peak_uniq$V2,"_",three_germ_peak_uniq$V3)
rownames(three_germ_peak_uniq) <- three_germ_peak_uniq$V5
head(three_germ_peak_uniq)

cluster_order <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/whole_Tss12000_v3/271_subcluster_matrix/cluster_ordered_meta_maxsample_new.txt")
colnames(ct_mt) <- cluster_order$germ_ratio_celltype

ectoderm_celltype <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/ectoderm_celltype.txt")  ## 91
endoderm_celltype <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/endoderm_celltype.txt")  ## 29
mesoderm_celltype <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/mesoderm_celltype.txt") ## 124

E7.5_ectoderm <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/E7.5_ectoderm_celltype.txt")
E7.5_endoderm <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/E7.5_endoderm_celltype.txt")
E7.5_mesoderm <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/E7.5_mesoderm_celltype.txt")
name <- c(as.character(E7.5_ectoderm$x),as.character(E7.5_mesoderm$x),as.character(E7.5_endoderm$x),as.character(ectoderm_celltype$x), as.character(mesoderm_celltype$x),as.character(endoderm_celltype$x))

peak_specific <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/sub_peak_ov_master_peak/celltype_peak_ov_whole/wholemap_germ_specific_OV_germlayer_celltype_add.bed")

remain_peak <- setdiff(as.character(three_germ_peak_uniq$V5),as.character(peak_specific$V5))
length(remain_peak)

germ_peak_all_cluster <- ct_mt[as.character(remain_peak), name]
peak_cluster_num <- rowSums(germ_peak_all_cluster)
peak_nocluster <- which(peak_cluster_num <= 3) 
peakname_nocluster <- rownames(germ_peak_all_cluster)[peak_nocluster]
germ_nocluster <- three_germ_peak_uniq[as.character(peakname_nocluster),]

ect_peak <- c(as.character(peak_specific$V5)[1:659], as.character(germ_nocluster$V5)[1:1312], setdiff(as.character(three_germ_peak_uniq$V5)[1:4715],c(as.character(peak_specific$V5)[1:659], as.character(germ_nocluster$V5)[1:1312])))
setdiff(ect_peak, as.character(three_germ_peak_uniq$V5)[1:4715])

mes_peak <- c(as.character(peak_specific$V5)[660:1147], as.character(germ_nocluster$V5)[1313:2098], setdiff(as.character(three_germ_peak_uniq$V5)[4716:8691],c(as.character(peak_specific$V5)[660:1147], as.character(germ_nocluster$V5)[1313:2098])))
setdiff(mes_peak, as.character(three_germ_peak_uniq$V5)[4716:8691])

end_peak <- c(as.character(peak_specific$V5)[1148:1353], as.character(germ_nocluster$V5)[2099:2521], setdiff(as.character(three_germ_peak_uniq$V5)[8692:10511],c(as.character(peak_specific$V5)[1148:1353], as.character(germ_nocluster$V5)[2099:2521])))
setdiff(end_peak, as.character(three_germ_peak_uniq$V5)[8692:10511])

### save peak info ###
### maintain peak ###
write.table(as.character(peak_specific$V5)[1:659], file = "Ectoderm_germlayer_maintain_peak.txt", quote = F, row.names = F, col.names = F)
write.table(as.character(peak_specific$V5)[660:1147], file = "Mesoderm_germlayer_maintain_peak.txt", quote = F, row.names = F, col.names = F)
write.table(as.character(peak_specific$V5)[1148:1353], file = "Endoderm_germlayer_maintain_peak.txt", quote = F, row.names = F, col.names = F)

### closed peak ###
write.table(as.character(germ_nocluster$V5)[1:1312], file = "Ectoderm_germlayer_closed_peak.txt", quote = F, row.names = F, col.names = F)
write.table(as.character(germ_nocluster$V5)[1313:2098], file = "Mesoderm_germlayer_closed_peak.txt", quote = F, row.names = F, col.names = F)
write.table(as.character(germ_nocluster$V5)[2099:2521], file = "Endoderm_germlayer_closed_peak.txt", quote = F, row.names = F, col.names = F)

### universal peak ###
write.table(setdiff(as.character(three_germ_peak_uniq$V5)[1:4715],c(as.character(peak_specific$V5)[1:659], as.character(germ_nocluster$V5)[1:1312])), file = "Ectoderm_germlayer_universal_peak.txt", quote = F, row.names = F, col.names = F)
write.table(setdiff(as.character(three_germ_peak_uniq$V5)[4716:8691],c(as.character(peak_specific$V5)[660:1147], as.character(germ_nocluster$V5)[1313:2098])), file = "Mesoderm_germlayer_universal_peak.txt", quote = F, row.names = F, col.names = F)
write.table(setdiff(as.character(three_germ_peak_uniq$V5)[8692:10511],c(as.character(peak_specific$V5)[1148:1353], as.character(germ_nocluster$V5)[2099:2521])), file = "Endoderm_germlayer_universal_peak.txt", quote = F, row.names = F, col.names = F)

peak_reorder <- c(ect_peak,mes_peak,end_peak)

germ_peak_all_cluster <- ct_mt[as.character(three_germ_peak_uniq$V5), name]
ct_mt_order <- germ_peak_all_cluster[peak_reorder,]

options(expression = 500000)
library(Matrix)
library(ComplexHeatmap)

#par(mar=c(10,10,10,10))

setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/11_subcluster_peakcalling_bed/germlayer_specific/version5_only_germcell/ordered_germ_specific_peak")
pdf("1_celltype_germpeak_all_celltype_new.pdf",width = 80,height = 80)
htmap = Heatmap(as.matrix(ct_mt_order), cluster_columns = F, cluster_rows = F, show_row_names = F, col = c("white", "black"), column_dend_height = unit(20, "cm"), column_dend_side = "bottom")
htmap
dev.off()

