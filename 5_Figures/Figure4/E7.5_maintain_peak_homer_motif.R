#### homer result - known motif ####
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/homer_motif")

ct.ls <- c("E7.5_Ectoderm_maintain_MotifOutput", "E7.5_Mesoderm_maintain_MotifOutput", "E7.5_Endoderm_maintain_MotifOutput")
motif.file.ls <- list()
motif.score.ls <- list()
high50.ls <- list()
top5.motif.ls <- list()

for (i in ct.ls){
  motif.file.ls[[i]] <- read.table(paste0(i, "/knownResults.txt"), sep = "\t", comment.char = "", header = T)
  motif.score.ls[[i]] <- data.frame(motif = motif.file.ls[[i]]$Motif.Name, Log.P.value = motif.file.ls[[i]]$Log.P.value)
  idx <- which((-motif.score.ls[[i]]$Log.P.value) > 50)
  high50.ls[[i]] <- motif.score.ls[[i]][idx,]
}

high50_motif_df <- Reduce(rbind, high50.ls)
high50_motif <- unique(as.character(high50_motif_df$motif))

motif_sub.ls <- list()
motif_sub_order.ls <- list()
motif_sub_order_p.ls <- list()

for (i in ct.ls){
  motif_sub.ls[[i]] <- motif.score.ls[[i]][which(motif.score.ls[[i]]$motif %in% high50_motif),]
  rownames(motif_sub.ls[[i]]) <- as.character(motif_sub.ls[[i]]$motif)
  motif_sub_order.ls[[i]] <- motif_sub.ls[[i]][high50_motif,]
  motif_sub_order_p.ls[[i]] <- motif_sub.ls[[i]][high50_motif, "Log.P.value"]
}

high50_motif_df <- Reduce(cbind, motif_sub_order_p.ls)
rownames(high50_motif_df) <- high50_motif
colnames(high50_motif_df) <- gsub("_MotifOutput", "", ct.ls)

high50_motif_df
#### order ####
#high50_motif_df_p_sel <- high50_motif_df[rownames(TF_gene_inuse),]

########### homer gene-motif correlation #######################
TF_gene <- readRDS("homer_motif_gene_correlarion.rds")
remian_TF <- setdiff(high50_motif, as.character(TF_gene$motif))
gene1 <- c("Pou5f1", "Dlx3", "Rfx2", "Gsc", "Rfx1", "Isl1", "Xbp1", "Sox9", "Crx", "Tbx6", "Klf14", "Tbx21", "Eomes", "Tbx21", "Foxa1", "Foxa2", "Foxf1", "Foxk1", "Foxp1", "Foxk2", "Foxo1")
TF_gene2 <- data.frame(motif = remian_TF, gene = gene1)
TF_gene3 <- rbind(TF_gene, TF_gene2)
saveRDS(TF_gene3, file = "homer_motif_gene_correlarion.rds")
write.table(TF_gene3, file = "homer_motif_gene_correlarion.txt", quote = F, row.names = F, sep = "\t")

rownames(TF_gene3) <- as.character(TF_gene3$motif)
TF_gene3_order <- TF_gene3[rownames(high50_motif_df),]

idx <- which(TF_gene3_order$gene %in% c("Foxa1","Foxa2", "Tbx21"))
TF_gene_inuse <- TF_gene3_order[setdiff(1:55, c(4,idx[1:3])),]

high50_motif_df_p_sel <- high50_motif_df[rownames(TF_gene_inuse),]
############# add gene expression matrix #########################
ct_gene <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/groupList/gene_271_celltypes_expression.rds")
ct_gene_sel <- ct_gene[as.character(TF_gene_inuse$gene), 1:10]

E7.5_ect <- apply(ct_gene_sel[,1:4], 1, mean)
E7.5_mes <- apply(ct_gene_sel[,5:8], 1, mean)
E7.5_end <- apply(ct_gene_sel[,9:10], 1, mean)
ct_gene_sel_ave <- cbind(E7.5_ect,E7.5_mes,E7.5_end)

all(rownames(high50_motif_df_p_sel) == as.character(TF_gene_inuse$motif))

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

ct_gene_sel_zscore <- .rowZscores(ct_gene_sel_ave, limit = T)
#which(rownames(ct_gene_sel_zscore) %in% "Slc22a16") #3
#ct_gene_sel_zscore_sel <- ct_gene_sel_zscore[setdiff(1:52,3),]

motif_mt <- data.frame(type = c(rep("E7.5_ect",51), rep("E7.5_mes",51),rep("E7.5_end",51)), motif = rep(as.character(TF_gene_inuse$motif),3), E7.5_expression = as.vector(ct_gene_sel_zscore_sel), log10.p = as.vector(high50_motif_df_p_sel))

id <- strsplit(as.character(motif_mt$motif), "/")
name <- unlist(lapply(seq(id), function(x){id[[x]][1]}))

motif_mt$motif_ab <- name

motif_mt$motif_ab <- factor(motif_mt$motif_ab, levels = c(rev(name[1:51])))
motif_mt$type <- factor(motif_mt$type, levels = c("E7.5_ect","E7.5_mes","E7.5_end"))


pdf("E7.5_TF_high50_motif_expression_plot.pdf", width = 5.2, height = 10)
p<- ggplot(motif_mt, aes(type, motif_ab, color = E7.5_expression)) +
  geom_point(aes(size=-log10.p)) +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  scale_color_gradient(limits = c(-2,2), low = "white", high = "red")
p
dev.off()

