library(DESeq2)
count <- read.csv('count.csv',header=T,row.names = 1)
count1  <- count[,-c(8:9)]
group_list <- as.factor(rep(c('ES','FOS','HSF','MYCN','ME','MYC_ME','MYCN_ME','DE'),c(2,2,2,2,2,2,2,2)))
colData <- data.frame(row.names=colnames(count1[,-1]), group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = count1[,-1], colData = colData, design = ~ group_list)
##fpkm
mcols(dds)$basepairs <- count1$Length
my_fpkm <- fpkm(dds)
head(my_fpkm)
write.csv(my_fpkm,"fpkm_n.csv")

library(edgeR)
##edgeR
count <- read.csv('count.csv',header = T,row.names = 1)
count  <- count[,c(1:3,12:13,18:19)]
dim(count)
head(count)
a <- rowSums(count[,-1])
count <- count[a>0,]
dim(count)
count <- count[,-1]
group_list <- as.factor(rep(c('ES','ME','DE'),c(2,2,2)))
head(count)
y<- DGEList(counts=count,group=group_list,genes=rownames(count))
dge <- calcNormFactors(y)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
dge <- estimateDisp(dge, design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
sig_cut=0.05

####define ES ME DE signatures
######
A1 <- makeContrasts(A=ES-ME,levels = design)
A2 <- makeContrasts(A=ES-DE,levels = design)
B1 <- makeContrasts(A=ME-ES,levels = design)
B2 <- makeContrasts(A=ME-DE,levels = design)
C1 <- makeContrasts(A=DE-ES,levels = design)
C2 <- makeContrasts(A=DE-ME,levels = design)

lrtA1 <- glmQLFTest(fit, contrast=A1)
lrtA2 <- glmQLFTest(fit, contrast=A2)
lrtB1 <- glmQLFTest(fit, contrast=B1)
lrtB2 <- glmQLFTest(fit, contrast=B2)
lrtC1 <- glmQLFTest(fit, contrast=C1)
lrtC2 <- glmQLFTest(fit, contrast=C2)

diffA1<-topTags(lrtA1,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffA2<-topTags(lrtA2,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffB1<-topTags(lrtB1,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffB2<-topTags(lrtB2,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffC1<-topTags(lrtC1,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffC2<-topTags(lrtC2,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table

diffA1 <- diffA1[diffA1$logFC>1.5,]
diffA2 <- diffA2[diffA2$logFC>1.5,]
diffB1 <- diffB1[diffB1$logFC>1.5,]
diffB2 <- diffB2[diffB2$logFC>1.5,]
diffC1 <- diffC1[diffC1$logFC>1.5,]
diffC2 <- diffC2[diffC2$logFC>1.5,]
dim(diffA1)
dim(diffA2)
dim(diffB1)
dim(diffB2)
dim(diffC1)
dim(diffC2)
C0 <- Reduce(intersect,list(diffA1$genes,diffA2$genes))
length(C0)#1753
C1 <- rbind(diffA1,diffA2)
dim(C1) # 5723    6
C1_s <- C1[C1$genes%in%C0,]
dim(C1_s)##3506    6
logfc <- aggregate(C1_s$logFC,by =list(C1_s$genes),min)
p_value <- aggregate(C1_s$PValue,by =list(C1_s$genes),max)
fdr <- aggregate(C1_s$FDR,by =list(C1_s$genes),max)
head(logfc)
head(fdr)
head(p_value)
sig_c0 <- cbind(logfc,p_value[,2],fdr[,2])
head(sig_c0)
colnames(sig_c0) <- c('gene_id','logfc','p_value','fdr')
dim(sig_c0)
sig_c0s <- sig_c0[order(sig_c0$logfc,decreasing = T),]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$logfc>1.5,]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$fdr<0.05,]
dim(sig_c0s)
write.csv(sig_c0s,'ES_marker.csv')

##ME marker
C0 <- Reduce(intersect,list(diffB1$genes,diffB2$genes))
length(C0)#1753
C1 <- rbind(diffB1,diffB2)
dim(C1) # 5723    6
C1_s <- C1[C1$genes%in%C0,]
dim(C1_s)##3506    6
logfc <- aggregate(C1_s$logFC,by =list(C1_s$genes),min)
p_value <- aggregate(C1_s$PValue,by =list(C1_s$genes),max)
fdr <- aggregate(C1_s$FDR,by =list(C1_s$genes),max)
head(logfc)
head(fdr)
head(p_value)
sig_c0 <- cbind(logfc,p_value[,2],fdr[,2])
head(sig_c0)
colnames(sig_c0) <- c('gene_id','logfc','p_value','fdr')
dim(sig_c0)
sig_c0s <- sig_c0[order(sig_c0$logfc,decreasing = T),]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$logfc>1.5,]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$fdr<0.05,]
dim(sig_c0s)
write.csv(sig_c0s,'ME_marker.csv')

##DE marker
C0 <- Reduce(intersect,list(diffC1$genes,diffC2$genes))
length(C0)
C1 <- rbind(diffC1,diffC2)
dim(C1) 
C1_s <- C1[C1$genes%in%C0,]
dim(C1_s)
logfc <- aggregate(C1_s$logFC,by =list(C1_s$genes),min)
p_value <- aggregate(C1_s$PValue,by =list(C1_s$genes),max)
fdr <- aggregate(C1_s$FDR,by =list(C1_s$genes),max)
head(logfc)
head(fdr)
head(p_value)
sig_c0 <- cbind(logfc,p_value[,2],fdr[,2])
head(sig_c0)
colnames(sig_c0) <- c('gene_id','logfc','p_value','fdr')
dim(sig_c0)
sig_c0s <- sig_c0[order(sig_c0$logfc,decreasing = T),]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$logfc>1.5,]
dim(sig_c0s)
sig_c0s <- sig_c0s[sig_c0s$fdr<0.05,]
dim(sig_c0s)
write.csv(sig_c0s,'DE_marker.csv')



rm(list = ls())
######
library(pheatmap)
data <- read.csv('fpkm_n.csv',header = T)
head(data)
dim(data)
colnames(data)[1] <- 'gene'
ES_d <- data[data$gene%in%c('POU5F1','SOX2','NANOG'),]
ME_d <- data[data$gene%in%c('MIXL1','GSC','EOMES'),]
dim(ME_d)
DE_d  <- data[data$gene%in%c('SOX17','FOXA2'),]
dim(DE_d)
marker_d <- rbind(ES_d,ME_d,DE_d)
head(marker_d)
marker_d$ES = rowMeans(marker_d[,2:3])
marker_d$FOS = rowMeans(marker_d[,4:5])
marker_d$HSF1 = rowMeans(marker_d[,6:7])
marker_d$MYCN = rowMeans(marker_d[,8:9])
marker_d$ME = rowMeans(marker_d[,10:11])
marker_d$MYC_ME = rowMeans(marker_d[,12:13])
marker_d$MYCN_ME = rowMeans(marker_d[,14:15])
marker_d$DE = rowMeans(marker_d[,16:17])
dim(marker_d)
head(marker_d)
marker_d <- marker_d[,c(1,18:24)]
head(marker_d)
library(ggplot2)
rownames(marker_d) <- marker_d$gene
marker_d <- marker_d[,-1]
head(marker_d)
class(marker_d)
d= t(marker_d)
class(d)
d <- as.data.frame(d)
d$group <- rownames(d)
d$group <- factor(d$group,levels = d$group)
head(d)
ggplot() + geom_bar(data = d, aes(x = group, y = POU5F1), 
                    stat = "identity")+theme_bw()


ES_marker <- read.csv('ES_marker.csv',header = T,row.names = 1)
ME_marker <- read.csv('ME_marker.csv',header = T,row.names = 1)
DE_marker <- read.csv('DE_marker.csv',header = T,row.names = 1)
gene_list <- data.frame(gene= c(ES_marker$gene_id,ME_marker$gene_id,DE_marker$gene_id))
head(gene_list)
data1 <- merge(gene_list,data,all=F,sort=F)
head(data1)
data2 <- data1[,c(1,2,3,10,11,16,17)]
head(data2)
rownames(data2)<- data2$gene
head(data2)
data2 <- data2[,-1]
annotation_col = data.frame(
  state = factor(rep(c("ES", "ME","DE"), each=2)))
rownames(annotation_col) = colnames(data2)
head(annotation_col)
ann_colors = list(
  state = c(ES = "mistyrose3", ME = "palegreen4", DE = "slategray1")
)
head(ann_colors)
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE, annotation_col=annotation_col,annotation_colors = ann_colors,labels_row = '')


#top 20
gene_list <- data.frame(gene= c(ES_marker$gene_id[1:20],ME_marker$gene_id[1:20],DE_marker$gene_id[1:20]))
data1 <- merge(gene_list,data,all=F,sort=F)
head(data1)
data2 <- data1[,c(1,2,3,10,11,16,17)]
head(data2)
rownames(data2)<- data2$gene
head(data2)
data2 <- data2[,-1]
annotation_col = data.frame(
  state = factor(rep(c("ES", "ME","DE"), each=2)))
rownames(annotation_col) = colnames(data2)
head(annotation_col)
ann_colors = list(
  state = c(ES = "mistyrose3", ME = "palegreen4", DE = "slategray1")
)
head(ann_colors)
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE, 
         annotation_col=annotation_col,annotation_colors = ann_colors,cellwidth = 50, cellheight = 10,main = 'Top20')


##top20 expression level across all samples
data <- read.csv('fpkm_n.csv',header = T)
colnames(data)[1] <- 'gene'
ES_marker <- read.csv('ES_marker.csv',header = T,row.names = 1)
ME_marker <- read.csv('ME_marker.csv',header = T,row.names = 1)
DE_marker <- read.csv('DE_marker.csv',header = T,row.names = 1)
gene_list <- data.frame(gene= c(ES_marker$gene_id[1:20],ME_marker$gene_id[1:20]))
gene_list1 <- gene_list
gene_list$gene[21] <- gene_list1$gene[28]
gene_list$gene[28] <- gene_list1$gene[21]
gene_list$gene[22] <- gene_list1$gene[35]
gene_list$gene[35] <- gene_list1$gene[22]
data1 <- merge(gene_list,data,all=F,sort=F)
head(data1[,1:4])
rownames(data1)<- data1$gene
head(data1[,1:6])
#data2 <- data1[,c(1,4:13)]
#head(data2)
#rownames(data2)<- data2$gene
#head(data2)
#data2 <- data2[,-1]
data2 <- data1[,-1]
data2 <- data2[,3:10]
head(data2)
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE,cellwidth = 20, cellheight = 10)


#
gene_list <- data.frame(gene= c(ME_marker$gene_id[1:20],DE_marker$gene_id[1:20]))
gene_list1 <- gene_list
gene_list1$gene[2:11] <- gene_list$gene[c(8,9,11:18)]
gene_list1$gene[12:13]<- gene_list$gene[3:4]
gene_list1$gene[14:16]<- gene_list$gene[c(19,2,20)]
gene_list1$gene[17:20]<- gene_list$gene[c(5,6,7,10)]
gene_list1$gene[21:24] <-  gene_list$gene[c(24,27,39:40)]
gene_list1$gene[25:40] <-  gene_list$gene[c(21,22,23,25,26,28:38)]
data1 <- merge(gene_list1,data,all=F,sort=F)
head(data1)
rownames(data1)<- data1$gene
head(data1)
data2 <- data1[,c(1,12:17)]
head(data2)
rownames(data2)<- data2$gene
head(data2)
data2 <- data2[,-1]
head(data2)
data3 <- data2[,-c(3:4)]
pheatmap(data2, scale = 'row', cluster_rows = FALSE,cluster_cols =   FALSE,cellwidth = 20, cellheight = 10)

##PCA and t-SNE
library(Rtsne)
library(ggplot2)
data <- read.csv('fpkm_n.csv',header = T,row.names = 1)
head(data)
data  <- data[,-c(11,12)]
pca <- prcomp(t(data))
pca_root <- pca$x
head(pca_root)
set.seed(1000000)
tsne_out <- Rtsne(
  pca_root[,1:14],
  dims = 2,
  pca = FALSE,
  perplexity = 4,
  max_iter = 10000
)
tsne_root <- tsne_out$Y
head(tsne_root)
group <- rep(c('ES','FOS','HSF1','MYCN','ME','MYCN_ME','DE'),c(2,2,2,2,2,2,2))
group <- factor(group,levels = c('ES','FOS','HSF1','MYCN','ME','MYCN_ME','DE'))
x = data.frame(tsne_root, group = group )
head(x)
colnames(x)[1:2] <- c('tSNE_1','tSNE_2')
ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(size = 8)  + theme_bw() 

y  = data.frame(pca_root[,1:2], group = group )
colnames(y)[1:2] <- c('PC1','PC2')
ggplot(y, aes(x=PC1, y=PC2, color=group)) + geom_point(size = 5)  + theme_bw() 

###correlation
library(pheatmap)
data <- read.csv('fpkm_n.csv',header = T,row.names =1 )
head(data)
dim(data)
data <- data[,-c(11,12)]
head(data)
data <- data[rowMeans(data)>1,]
dim(data)
c <- cor(data)
dim(c)
pheatmap(c)


library(ggplot2)
data <- read.csv('fpkm_n.csv',header = T,row.names = 1)
head(data)
OCT4 <- data[rownames(data)=='OCT4',] 
SOX2 <- data[rownames(data)=='SOX2',] 
MIXL1 <- data[rownames(data)=='MIXL1',]
group <- rep(c('ES','FOS','HSF1','MYCN','ME','MYC_ME','DE'),c(2,2,2,2,2,2,2))
group <- factor(group,levels = c('ES','FOS','HSF1','MYCN','ME','MYC_ME','DE')) 
exp <- data.frame(SOX2= as.numeric(SOX2),MIXL1 =as.numeric(MIXL1),group=group)
head(exp)
ggplot(exp, aes(x=group, y=SOX2, color = group)) +
  geom_boxplot()
ggplot(exp, aes(x=group, y=MIXL1, color = group)) +
  geom_boxplot()
