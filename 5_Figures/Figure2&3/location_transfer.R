library(Seurat)
library(Matrix)
##CCA co-embedding
data1=readRDS("Geo.rds")
data2=readRDS("genescore_mat_seurat.rds")
mat1=data1@assays$RNA@counts
mat2=data2@assays$RNA@counts
rowmeans1=rowMeans(mat1)
rowmeans2=rowMeans(mat2)
idx1=which(rowmeans1 != 0)
idx2=which(rowmeans2 != 0)
mat1=mat1[idx1,]
mat2=mat2[idx2,]
gene_name1=rownames(mat1)
gene_name2=rownames(mat2)
gene_intersect=intersect(gene_name1,gene_name2)
data_cca <- RunCCA(object1 = data1, object2 = data2,features=gene_intersect,num.cc = 20)
saveRDS(data_cca,file="after_cca.rds")

##L2 normalization
data_L2=L2Dim(object=data_cca,reduction="cca")

##calculate Euclid distance between co-embedding ATAC cells and locations
E_dist=function(a,b)
{
    p=(a-b)
    p=sum(p*p)
    p=sqrt(p)
    return(p)
}


mat=data_L2@reductions$cca.l2@cell.embeddings
nrow=nrow(mat)
ncol=nrow-83
similar=matrix(data=NA, nrow = 83, ncol = ncol, byrow = FALSE, dimnames = NULL)

for (i in 1:83)
{   
    for (j in 84:nrow)
    {
     	similar[i,j-83]=E_dist(mat[i,],mat[j,])
    }
}

##select clostest position to ATAC cell
similar=t(similar)
nrow=nrow(similar)
ncol=ncol(similar)
select=matrix(data=0, nrow = nrow, ncol = ncol, byrow = FALSE, dimnames = NULL)
threshold=1 ##threshold means maxium distance
for (i in 1:nrow)
{
    gene=similar[i,]
    if(min(gene)<threshold)
    {
     	j=which(gene==min(gene))
        select[i,j]=1
    }
    
}
CCA=readRDS("after_cca.rds")
mat=CCA@assays$RNA@counts
name=colnames(mat)
rowname=name[84:length(name)]
colname=name[1:83]
rownames(select)=rowname
colnames(select)=colname
saveRDS(select,file="select_L2.rds")
write.table(select,file="select_L2.txt")

##get genescore V.S.location  matrix
sc_mat=data2@assays$RNA@data
bu_mat=read.table("select_L2.txt")
bu_mat=as.matrix(bu_mat)
inte_mat=sc_mat %*% bu_mat
w=colSums(bu_mat)
w=as.numeric(w)
for (i in 1:ncol(bu_mat))
{
    inte_mat[,i]=inte_mat[,i]/w[i]
}
saveRDS(inte_mat,file="inte_mat_genescore.rds")

##get accessibility V.S.location  matrix
sc_mat=readRDS("E7_5_germlayer_peak.rds")
bu_mat=read.table("select_L2.txt")
bu_mat=as.matrix(bu_mat)
inter_mat=sc_mat %*% bu_mat
w=colSums(bu_mat)
w=as.numeric(w)
for (i in 1:ncol(bu_mat))
{
    inter_mat[,i]=inter_mat[,i]/w[i]
}
saveRDS(inter_mat,file="inter_mat_ATAC.rds")

