library(Rtsne)
library(ggplot2)
data <- read.csv('fpkm.csv',header = T,row.names = 1)
pca <- prcomp(t(data))
pca_root <- pca$x
head(pca)
set.seed(1000000)
tsne_out <- Rtsne(
  pca_root[,1:13],
  dims = 2,
  pca = FALSE,
  perplexity = 4,
  max_iter = 10000
)
tsne_root <- tsne_out$Y
head(tsne_root)
group <- rep(c('ES','FOS','HSF','MYCN','ME','MYC_ME','DE'),c(2,2,2,2,2,2,2))
group <- factor(group,levels = c('ES','FOS','HSF','MYCN','ME','MYC_ME','DE'))
x = data.frame(tsne_root, group = group )
head(x)
colnames(x)[1:2] <- c('tSNE_1','tSNE_2')
ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(size = 8)  + theme_bw() 

y  = data.frame(pca_root[,1:2], group = group )
colnames(y)[1:2] <- c('PC1','PC2')
ggplot(y, aes(x=PC1, y=PC2, color=group)) + geom_point(size = 5)  + theme_bw() 
rm(list = ls())
###
library(edgeR)
##edgeR
data <- read.csv('fpkm.csv',header = T,row.names = 1)
dim(data)
head(data)
##ES
ES_fc1 <- rowMeans(data[,1:2])/rowMeans(data[,9:10]) 
ES_fc2 <- rowMeans(data[,1:2])/rowMeans(data[,13:14]) 
ES_fc <- data.frame(gene= names(ES_fc1),fc1=ES_fc1,fc2=ES_fc2)
head(ES_fc)
fc <- apply(ES_fc[,2:3],1,min)
fc
ES_fc$fc <- fc
ES_fc = ES_fc[order(ES_fc$fc,decreasing = T),]
head(ES_fc)
ES_fc$logfc <- log2(ES_fc$fc)
write.csv(ES_fc,'ES_sig.csv')

##ME
ME_fc1 <- rowMeans(data[,9:10])/rowMeans(data[,1:2]) 
ME_fc2 <- rowMeans(data[,9:10])/rowMeans(data[,13:14]) 
ME_fc <- data.frame(gene= names(ME_fc1),fc1=ME_fc1,fc2=ME_fc2)
head(ME_fc)
fc <- apply(ME_fc[,2:3],1,min)
fc
ME_fc$fc <- fc
ME_fc$logfc <- log2(ME_fc$fc)
write.csv(ME_fc,'ME_sig.csv')

##DE
DE_fc1 <- rowMeans(data[,13:14])/rowMeans(data[,1:2]) 
DE_fc2 <- rowMeans(data[,13:14])/rowMeans(data[,9:10]) 
DE_fc <- data.frame(gene= names(DE_fc1),fc1=DE_fc1,fc2=DE_fc2)
head(DE_fc)
fc <- apply(DE_fc[,2:3],1,min)
fc
DE_fc$fc <- fc
DE_fc$logfc <- log2(DE_fc$fc)
write.csv(DE_fc,'DE_sig.csv')

rm(list = ls())
######
library(pheatmap)
data <- read.csv('fpkm.csv',header = T)
head(data)
colnames(data)[1] <- 'gene'
dim(data)
ES_marker <- read.csv('ES_sig.csv',header = T,row.names = 1)
ME_marker <- read.csv('ME_sig.csv',header = T,row.names = 1)
DE_marker <- read.csv('DE_sig.csv',header = T,row.names = 1)

##ES ME DE
gene_list <- data.frame(gene= c(ES_marker$gene[1:20],ME_marker$gene[1:20],DE_marker$gene[1:20]))
gene_list1 <- gene_list
data1 <- merge(gene_list,data,all=F,sort=F)
head(data1[,1:4])
dim(data1)
rownames(data1)<- data1$gene
head(data1)
data2 <- data1[,c(2:3,10:11,14:15)]
head(data2)
#rownames(data2)<- data2$gene
#head(data2)
#data2 <- data2[,-1]
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE)


##ES ME
a <- ME_marker$gene[1:20]
a1 <- a[c(15,16,20)]
a2 <- a[-c(15,16,20)]
a = c(a1,a2)
gene_list <- data.frame(gene= c(ES_marker$gene[1:20],a))
gene_list1 <- gene_list
data1 <- merge(gene_list,data,all=F,sort=F)
head(data1[,1:4])
rownames(data1)<- data1$gene
head(data1)
data2 <- data1[,c(4:11)]
head(data2)
#rownames(data2)<- data2$gene
#head(data2)
#data2 <- data2[,-1]
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE)

## ME  DE
gene_list <- data.frame(gene= c(ME_marker$gene[1:20],DE_marker$gene[1:20]))
gene_list1 <- gene_list
gene_list1[21,1] <- gene_list[25,1]
gene_list1[22,1] <- gene_list[34,1]
gene_list1[23,1] <- gene_list[36,1]
gene_list1[24,1] <- gene_list[22,1]
gene_list1[34,1] <- gene_list[24,1]
gene_list1[36,1] <- gene_list[23,1]


data1 <- merge(gene_list1,data,all=F,sort=F)
head(data1[,1:4])
rownames(data1)<- data1$gene
head(data1)
data2 <- data1[,c(12:15)]
head(data2)
#rownames(data2)<- data2$gene
#head(data2)
#data2 <- data2[,-1]
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =  FALSE)
##pearson
###correlation
library(pheatmap)
data <- read.csv('fpkm.csv',header = T,row.names =1 )
head(data)
dim(data)
data <- data[rowMeans(data)>1,]
dim(data)
c <- cor(data)
dim(c)
pheatmap(c)


