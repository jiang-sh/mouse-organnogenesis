TF_gene_cor <- readRDS("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/TF_expression/TF_gene_cor.rds")

############################ 1. motif enrichment #################################################
setwd("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/17_TF_analysis/fetal_celltype_peak_motif/")

ct.motif <- read.table("motif_file_ordered.txt")
ct.ls <- as.character(ct.motif$V1)

motif.file.ls <- list()
motif.score.ls <- list()

for (i in ct.ls){
  motif.file.ls[[i]] <- read.table(paste0(i, "/knownResults.txt"), sep = "\t", comment.char = "", header = T)
  motif.score.ls[[i]] <- data.frame(motif = motif.file.ls[[i]]$Motif.Name, Log.P.value = motif.file.ls[[i]]$Log.P.value)
  #idx <- which((-motif.score.ls[[i]]$Log.P.value) > 100)
  #high100.ls[[i]] <- motif.score.ls[[i]][idx,]
  #top5.high00.ls[[i]] <- high100.ls[[i]][1:5,]
}

motif.order.ls <- list()
motif_order <- unique(as.character(motif.score.ls[[1]]$motif))

for (i in ct.ls){
  motif.file.ls[[i]] <- read.table(paste0(i, "/knownResults.txt"), sep = "\t", comment.char = "", header = T)
  motif <- data.frame(motif = motif.file.ls[[i]]$Motif.Name, Log.P.value = motif.file.ls[[i]]$Log.P.value)
  idx <- which(motif$motif %in% 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer')
  idy <- which(motif$motif %in% 'RORgt(NR)/EL4-RORgt.Flag-ChIP-Seq(GSE56019)/Homer')
  motif_sub <- motif[setdiff(1:428, c(idx[2], idy[2])),]

  rownames(motif_sub) <- as.character(motif_sub$motif)
  motif.order.ls[[i]] <- motif_sub[motif_order,]
}

p_value_ls <- lapply(motif.order.ls, function(x){x$Log.P.value})
ct_enrich <- Reduce(cbind, p_value_ls)
rownames(ct_enrich) <- rownames(motif.order.ls[[1]])
colnames(ct_enrich) <- names(p_value_ls)
ct_enrich_score <- -ct_enrich

saveRDS(ct_enrich_score, file = "TF_ct_428.rds")

motifname <- rownames(ct_enrich_score)
id <- strsplit(as.character(motifname), "/")
name <- unlist(lapply(seq(id), function(x){id[[x]][1]}))
motifname_ab <- name
rep_motif <- names(table(name)[which(table(name)>1)])
rep_score <- ct_enrich_score[sort(motifname[which(name %in% rep_motif)]),]
apply(rep_score, 1, summary)
rm_motif <- c("COUP-TFII(NR)/Artia-Nr2f2-ChIP-Seq(GSE46497)/Homer", "FOXA1(Forkhead)/LNCAP-FOXA1-ChIP-Seq(GSE27824)/Homer","GRE(NR),IR3/A549-GR-ChIP-Seq(GSE32465)/Homer","HRE(HSF)/Striatum-HSF1-ChIP-Seq(GSE38000)/Homer", "Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer","OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer","STAT6(Stat)/CD4-Stat6-ChIP-Seq(GSE22104)/Homer","THRb(NR)/HepG2-THRb.Flag-ChIP-Seq(Encode)/Homer","c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer","p53(p53)/Saos-p53-ChIP-Seq(GSE15780)/Homer","p53(p53)/mES-cMyc-ChIP-Seq(GSE11431)/Homer")

ct_enrich_score_sel <- ct_enrich_score[setdiff(1:426, which(rownames(ct_enrich_score) %in% rm_motif)),]

motifname2 <- rownames(ct_enrich_score_sel)
id <- strsplit(as.character(motifname2), "/")
name2 <- unlist(lapply(seq(id), function(x){id[[x]][1]}))
motifname_ab2 <- name2
table(table(motifname_ab2))

ov_motif <- intersect(TF_gene_cor$TF_ab, motifname_ab2)
gene1 <- TF_gene_cor[which(TF_gene_cor$TF_ab %in% ov_motif),] ## 304 

#not_ov_motif <- setdiff(motifname_ab2, TF_gene_cor$TF_ab)

##################### 2. gene celltype expression ################################
load("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/00_yuhao_peak2gene_link/gene_expression_normalize/1_average_expression.RData")
y <- t(as.matrix(cluster.average@assays$RNA@data))

intersect_gene <- intersect(colnames(y), gene1$gene)
gene1[which(gene1$gene %in% setdiff(gene1$gene, colnames(y))),]
gene1[120,2] <- "Nr4a1" 

intersect_gene <- intersect(colnames(y), gene1$gene) ## 291 gene
motif_celltype <- TF_gene_cor[which(TF_gene_cor$gene %in% intersect_gene),]  ## 303 motif 
rep_gene <- names(table(motif_celltype$gene)[table(motif_celltype$gene) > 1])
rep_gene_motif <- motif_celltype[which(motif_celltype$gene %in% rep_gene),]
rep_gene_motif_sort <- rep_gene_motif[order(rep_gene_motif$gene),]
rep_enrich <-  ct_enrich_score[which(rownames(ct_enrich_score) %in% rep_gene_motif$TF),]
rm_gene <- c("Ets1-distal(ETS)/CD4+-PolII-ChIP-Seq(Barski_et_al.)/Homer", "Unknown-ESC-element(?)/mES-Nanog-ChIP-Seq(GSE11724)/Homer","OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer","OCT:OCT(POU,Homeobox)/NPC-Brn1-ChIP-Seq(GSE35496)/Homer","OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)/mES-Oct4-ChIP-Seq(GSE11431)/Homer","RFX(HTH)/K562-RFX3-ChIP-Seq(SRA012198)/Homer", "Oct4:Sox17(POU,Homeobox,HMG)/F9-Sox17-ChIP-Seq(GSE44553)/Homer","Tbet(T-box)/CD8-Tbet-ChIP-Seq(GSE33802)/Homer","Fox:Ebox(Forkhead,bHLH)/Panc1-Foxa2-ChIP-Seq(GSE47459)/Homer", "GATA(Zf),IR3/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer", "FOXA1(Forkhead)/LNCAP-FOXA1-ChIP-Seq(GSE27824)/Homer","NF1:FOXA1(CTF,Forkhead)/LNCAP-FOXA1-ChIP-Seq(GSE27824)/Homer")
check_tf <- c("FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(GSE26831)/Homer", "FOXA1(Forkhead)/LNCAP-FOXA1-ChIP-Seq(GSE27824)/Homer","NF1:FOXA1(CTF,Forkhead)/LNCAP-FOXA1-ChIP-Seq(GSE27824)/Homer")
apply(ct_enrich_score[check_tf,],1, summary)

remain_tf <- setdiff(as.character(motif_celltype$TF), rm_gene) 
motif_celltype_sel <- motif_celltype[which(as.character(motif_celltype$TF) %in% remain_tf),]
celltype_enrich <- ct_enrich_score[which(rownames(ct_enrich_score) %in% motif_celltype_sel$TF),]
rename <- gsub("_MotifOutput", "", colnames(celltype_enrich))
colnames(celltype_enrich) <- rename
saveRDS(celltype_enrich, "selected_291_TFs_celltype_enrichment.rds")

celltype_gene <- y[,which(colnames(y) %in% motif_celltype_sel$gene)]
label_new <- read.table("/xtdisk/jiangl_group/yuchw/Mouse_Organogenesis_v2/whole_Tss12000_v3/271_subcluster_matrix/celltype_germ_maxratio_meta.txt")
label_rename <- unique(label_new)
label_rename_order <- label_rename[rownames(celltype_gene),]
rownames(celltype_gene) <- label_rename_order$germ_ratio_celltype

rownames(motif_celltype_sel) <- motif_celltype_sel$TF
motif_celltype_sel_order <- motif_celltype_sel[rownames(celltype_enrich),]
saveRDS(motif_celltype_sel_order, file = "selected_291_TF_genes_correlation.rds")

celltype_gene_order <- celltype_gene[colnames(celltype_enrich), motif_celltype_sel_order$gene]
celltype_expression <- t(celltype_gene_order)
saveRDS(celltype_expression, file = "selected_291_genes_celltype_expression.rds")
save.image("co-TF_finding.RData")


ct.exp.enrich <- list()
for (i in 1:length(rownames(celltype_expression))){
  exp <- as.vector(celltype_expression[i,])
  enrich <- as.vector(celltype_enrich[i,])
  ct.exp <- colnames(celltype_expression)[which(exp > 1 & exp > 2*mean(exp))]
  ct.tf <- colnames(celltype_enrich)[which(enrich > 100)]
  ct.exp.enrich[[rownames(celltype_expression)[i]]] <- intersect(ct.exp,ct.tf)
}

ct.exp.enrich.sel <- ct.exp.enrich[lapply(ct.exp.enrich, function(x){length(x)})>0]   ## 97
saveRDS(ct.exp.enrich.sel, file = "TF_expression_and_enrich_in_celltype.rds")

gene_name <- names(ct.exp.enrich.sel)
rownames(motif_celltype_sel_order) <- motif_celltype_sel_order$gene
motif_celltype_sel_order2 <- motif_celltype_sel_order[gene_name,]

tf_name <- as.character(motif_celltype_sel_order2$TF)

### 2. correlation ###
celltype_enrich_exp <- t(celltype_enrich[tf_name,])

library(Matrix)
ct2motif_scale = Matrix(0, ncol = 97, nrow = 271, sparse = T)
rownames(ct2motif_scale) <- rownames(celltype_enrich_exp)
colnames(ct2motif_scale) <- colnames(celltype_enrich_exp)

for (i in colnames(celltype_enrich_exp)){
  mt = celltype_enrich_exp[,i]
  ct2motif_scale[,i] <- (mt-min(mt))/(max(mt)-min(mt))*100
}

cor_mat = cor(as.matrix(ct2motif_scale), method="spearman")  ## 97, 97 

pdf("2_expression_TF_correlation_pheatmap.pdf", height = 10, width = 10)
p <- pheatmap::pheatmap(cor_mat[1:97,1:97]/apply(cor_mat[1:97,1:97], 1, max), fontsize=5)
p
dev.off()

