#### homer result - known motif ####
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/fetal_celltype_peak_motif/")

ct.motif <- read.table("motif_file_ordered.txt")
ct.ls <- as.character(ct.motif$V1)

########################## 3. plot maintain TF and their expression ##############
main_TF_exp <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/TF_expression/maintain_E7.5_motif_254_celltypes_expression.rds")
E13_mian_TF_and_germ_exp <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/TF_expression/TF_and_gene_E13.5_ect_mes_end.rds")

main_TF_ab <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/TF_expression/motif_name_order.rds")
rownames(E13_mian_TF_and_germ_exp) <- E13_mian_TF_and_germ_exp$motif_ab
E13_mian_TF_and_germ_exp_order <- E13_mian_TF_and_germ_exp[main_TF_ab,]
main_TF <- as.character(E13_mian_TF_and_germ_exp_order$motif)

motif_sub_main.ls <- list()
motif_sub_main_order.ls <- list()
motif_sub_main_order_p.ls <- list()

for (i in ct.ls){
  motif_sub_main.ls[[i]] <- motif.score.ls[[i]][which(motif.score.ls[[i]]$motif %in% main_TF),]
  rownames(motif_sub_main.ls[[i]]) <- as.character(motif_sub_main.ls[[i]]$motif)
  motif_sub_main_order.ls[[i]] <- motif_sub_main.ls[[i]][main_TF,]
  motif_sub_main_order_p.ls[[i]] <- motif_sub_main.ls[[i]][main_TF, "Log.P.value"]
}

main_motif_sum_df <- Reduce(cbind, motif_sub_main_order_p.ls)
colnames(main_motif_sum_df) <- gsub("_MotifOutput", "", ct.ls)

main_motif_sum_df_sel <- main_motif_sum_df[,1:244]
rownames(main_motif_sum_df_sel) <- main_TF
main_TF_exp_sel <- main_TF_exp[,11:254]
main_TF_exp_sel_order <- main_TF_exp[,colnames(main_motif_sum_df_sel)]
all(colnames(main_motif_sum_df_sel) == colnames(main_TF_exp_sel_order))

###################### 3.1 plot point TF and expression ################
celltype.ls <- unlist(lapply(colnames(main_motif_sum_df_sel),function(i){rep(i,34)}))
score <- as.vector(unlist(motif_sub_main_order_p.ls))[1:8296]
main_motif_df <- data.frame(motif = rep(main_TF, 244), type = celltype.ls, log10.p = score)

id <- strsplit(as.character(main_motif_df$motif), "/")
name <- unlist(lapply(seq(id), function(x){id[[x]][1]}))
main_motif_df$motif_ab <- name

main_motif_df$motif_ab <- factor(main_motif_df$motif_ab, levels = c(rev(main_TF_ab)))
main_motif_df$type <- factor(main_motif_df$type, levels = colnames(main_motif_sum_df_sel))

saveRDS(main_motif_df, file = "main_motif_TF_244_celltypes.rds")

main_TF_exp_row_zscore <- .rowZscores(main_TF_exp_sel_order, limit = T)
all(colnames(main_motif_sum_df_sel) == colnames(main_TF_exp_row_zscore))
main_ct_exp <- as.vector(main_TF_exp_row_zscore)
main_motif_df$expression <- main_ct_exp

########## bubble plot ##########3
library(ggplot2)
pdf("main_motif_and_TF_expression_plot_recolor_top5.pdf", width = 45, height = 7)
p <- ggplot(main_motif_df, aes(type, motif_ab, color = expression)) +
  geom_point(aes(size=-log10.p)) +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  scale_color_gradient2(limits = c(-2,2), low = "#4951AA", mid = "white", high = "#DD3423")
p
dev.off()

save.image("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/TF_expression/E13.5_celltypes_motif_and_expression.rds")

main_motif_df_sel <- main_motif_df[which(main_motif_df$motif_ab %in% TF_sel),]
main_motif_df_sel$motif <- factor(main_motif_df_sel$motif, levels = as.character(main_motif_df_sel$motif[1:15]))
main_motif_df_sel$motif_ab <- factor(main_motif_df_sel$motif_ab, levels = rev(as.character(main_motif_df_sel$motif_ab[1:15])))

#### learned select celltype ###
#main_motif_sum_df_sel2 <- main_motif_sum_df_sel[levels(main_motif_df_sel2$motif),]
main_motif_df_sel_new <- main_motif_df_sel
main_motif_df_sel_new$motif <- as.character(main_motif_df_sel_new$motif)
main_motif_df_sel_new$type <- as.character(main_motif_df_sel_new$type)
main_motif_df_sel_new$motif_ab <- as.character(main_motif_df_sel_new$motif_ab)

ect_sel <- main_motif_df_sel_new[which(main_motif_df_sel_new$motif %in% as.character(main_motif_df_sel_new$motif[1:5]) & main_motif_df_sel_new$type %in% type[1:91]),]
mes_sel <- main_motif_df_sel_new[which(main_motif_df_sel_new$motif %in% as.character(main_motif_df_sel_new$motif[6:10]) & main_motif_df_sel_new$type %in% type[92:215]),]
end_sel <- main_motif_df_sel_new[which(main_motif_df_sel_new$motif %in% as.character(main_motif_df_sel_new$motif[11:15]) & main_motif_df_sel_new$type %in% type[216:244]),]

type <- c("ect_sel","mes_sel","end_sel")
order <- list()
enrich <- list()
exp_enrich <- list()
for (i in type){
  order[[i]] <- sort(tapply(get(i)$expression, get(i)$type, median))
  enrich[[i]] <- sort(tapply(get(i)$log10.p, get(i)$type, median))
  tmp <- data.frame(type = names(order[[i]]), median_exp = as.vector(order[[i]]))
  tmp2 <- data.frame(type = names(enrich[[i]]), median_enrich = as.vector(enrich[[i]]))
  tmp3 <- merge(tmp, tmp2, by = "type")
  exp_enrich[[i]] <- tmp3[order(-tmp3$median_exp,tmp3$median_enrich), ]
}

ect_sel_exp_enrich <- exp_enrich[["ect_sel"]][which(exp_enrich[["ect_sel"]]$median_exp > 0 & exp_enrich[["ect_sel"]]$median_enrich < -500),][1:20,]
dim(ect_sel_exp_enrich)
mes_sel_exp_enrich <- exp_enrich[["mes_sel"]][which(exp_enrich[["mes_sel"]]$median_enrich < -1000),][1:20,]
dim(mes_sel_exp_enrich)
end_sel_exp_enrich <- exp_enrich[["end_sel"]][which(exp_enrich[["end_sel"]]$median_enrich < -120),]
dim(end_sel_exp_enrich)

celltype <- c(as.character(ect_sel_exp_enrich$type), as.character(mes_sel_exp_enrich$type), as.character(end_sel_exp_enrich$type))

main_motif_df_sel3 <- main_motif_df_sel[which(main_motif_df_sel$type %in% celltype),]
main_motif_df_sel3$type <- factor(main_motif_df_sel3$type, levels = celltype)

library(ggplot2)
pdf("main_motif_and_TF_expression_plot_recolor_select_high60ct_top5.pdf", width = 12, height = 4)
p <- ggplot(main_motif_df_sel3, aes(type, motif_ab, color = expression)) +
  geom_point(aes(size=-log10.p)) +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  scale_color_gradient2(limits = c(-2,2), low = "#4951AA", mid = "white", high = "#DD3423")
p
dev.off()



