##peak
##PCA analysis.
library(ggplot2)
count <- read.table('DNB_read_count_0.8.txt',header = F)
head(count)
sample <- rep(c('DE','ES','FOS','HSF','ME','MYC_ME','MYCN','TP53'),each=3)
s <- rep(c(1,2,3),9)
sample_name <-paste(sample,s,sep = '_')
sample_name
colnames(count)[5:31] <- sample_name
a <- paste(count$V1,count$V2,count$V3,sep = '_')
head(a)
rownames(count) <- a
data <- count[,-c(1:4)]
head(data)
dim(data)
data <- data[,-c(25:27)]
data_log <- log(data+1)
pca <- princomp(data_log)
pca_root <- pca$loadings[,1:2]
group <- as.factor(rep(c('DE','ES','FOS','HSF','ME','MYC_ME','MYCN','MYCN_ME'),each=3))
group
pca_data <- data.frame(pca_root,group = group)
head(pca_data)
##all samples
ggplot(data=pca_data, aes(x=Comp.1, y=Comp.2, color=group))+geom_point(size=5)+theme_bw()+xlab('PC1(24.06%)')+ylab('PC2(10.55%)')
##edgeR
library(edgeR)
group1 <-  as.factor(rep(c('DE','ES','ME'),each=3))
data1 <- data[,c(1:6,13:15)] 
head(data1)
y<-DGEList(counts=data1,group=group1,genes=rownames(data1))
dge <- calcNormFactors(y)
design <- model.matrix(~0+group1)
colnames(design) <- levels(group1)
dge <- estimateDisp(dge, design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
sig_cut=0.05
##ES
AB <- makeContrasts(A=ES-ME,levels = design)
AC <- makeContrasts(A=ES-DE,levels = design)
lrtAB <- glmQLFTest(fit, contrast=AB)
lrtAC <- glmQLFTest(fit, contrast=AC)
diffAB <-topTags(lrtAB,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffAC <-topTags(lrtAC,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
p <- intersect(diffAB$genes,diffAC$genes)
diffA <- rbind(diffAB[diffAB$genes%in%p,],diffAC[diffAC$genes%in%p,])
dim(diffA)
FDR <- aggregate(diffA$FDR,by=list(diffA$genes),max)
Pvalue <- aggregate(diffA$PValue,by=list(diffA$genes),max)
LogFC <- aggregate(diffA$logFC,by=list(diffA$genes),min)
ES <- cbind(FDR,Pvalue[,2],LogFC[,2])
head(ES)
colnames(ES) <- c('genes','FDR','Pvalue','LogFC')
head(ES)
ES <- ES[ES$FDR<=0.05,]
dim(ES)
ES <- ES[ES$LogFC>0.5,]
dim(ES)
write.csv(ES,'ES_sig_peaks.csv')


##ME
BA <- makeContrasts(A=ME-ES,levels = design)
BC <- makeContrasts(A=ME-DE,levels = design)
lrtBA <- glmQLFTest(fit, contrast=BA)
lrtBC <- glmQLFTest(fit, contrast=BC)
diffBA <-topTags(lrtBA,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffBC <-topTags(lrtBC,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
p <- intersect(diffBA$genes,diffBC$genes)
diffB <- rbind(diffBA[diffBA$genes%in%p,],diffBC[diffBC$genes%in%p,])
dim(diffB)
FDR <- aggregate(diffB$FDR,by=list(diffB$genes),max)
Pvalue <- aggregate(diffB$PValue,by=list(diffB$genes),max)
LogFC <- aggregate(diffB$logFC,by=list(diffB$genes),min)
ME <- cbind(FDR,Pvalue[,2],LogFC[,2])
head(ME)
colnames(ME) <- c('genes','FDR','Pvalue','LogFC')
head(ME)
ME <- ME[ME$FDR<=0.05,]
dim(ME)
ME <- ME[ME$LogFC>0.5,]
dim(ME)
write.csv(ME,'ME_sig_peaks.csv')
###DE
CA <- makeContrasts(A=DE-ES,levels = design)
CB <- makeContrasts(A=DE-ME,levels = design)
lrtCA <- glmQLFTest(fit, contrast=CA)
lrtCB <- glmQLFTest(fit, contrast=CB)
diffCA <-topTags(lrtCA,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffCB <-topTags(lrtCB,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
p <- intersect(diffCA$genes,diffCB$genes)
diffC <- rbind(diffCA[diffCA$genes%in%p,],diffCB[diffCB$genes%in%p,])
dim(diffC)
FDR <- aggregate(diffC$FDR,by=list(diffC$genes),max)
Pvalue <- aggregate(diffC$PValue,by=list(diffC$genes),max)
LogFC <- aggregate(diffC$logFC,by=list(diffC$genes),min)
DE <- cbind(FDR,Pvalue[,2],LogFC[,2])
head(DE)
colnames(DE) <- c('genes','FDR','Pvalue','LogFC')
head(DE)
DE <- DE[DE$FDR<=0.05,]
dim(DE)
DE <- DE[DE$LogFC>0.5,]
dim(DE)
write.csv(DE,'DE_sig_peaks.csv')

##Factor
library(edgeR)
library(ggplot2)
count <- read.table('DNB_read_count_0.8.txt',header = F)
head(count)
sample <- rep(c('DE','ES','FOS','HSF','ME','MYC_ME','MYCN','MYCN_ME','TP53'),each=3)
s <- rep(c(1,2,3),9)
sample_name <-paste(sample,s,sep = '_')
sample_name
colnames(count)[5:31] <- sample_name
a <- paste(count$V1,count$V2,count$V3,sep = '_')
head(a)
rownames(count) <- a
data <- count[,-c(1:4)]
head(data)
dim(data)
group <- as.factor(rep(c('DE','ES','FOS','HSF','ME','MYC_ME','MYCN','MYCN_ME','TP53'),each=3))
group
y<- DGEList(counts=data,group=group,genes=rownames(data))
dge <- calcNormFactors(y)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
sig_cut=0.01
B1 <- makeContrasts(A=FOS-ME,levels = design)
B2 <- makeContrasts(A=HSF-ME,levels = design)
B3 <- makeContrasts(A=TP53-ME,levels = design)
C1 <- makeContrasts(A=MYC_ME-DE,levels = design)
lrtB1 <- glmQLFTest(fit, contrast=B1)
lrtB2 <- glmQLFTest(fit, contrast=B2)
lrtB3 <- glmQLFTest(fit, contrast=B3)
lrtC1 <- glmQLFTest(fit, contrast=C1)

diffB1<-topTags(lrtB1,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffB2<-topTags(lrtB2,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffB3<-topTags(lrtB3,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffC1<-topTags(lrtC1,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffC2<-topTags(lrtC2,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
diffB1 <- diffB1[diffB1$logFC>2,]
diffB2 <- diffB2[diffB2$logFC>2,]
diffB3 <- diffB3[diffB3$logFC>2,]

diffC1 <- diffC1[diffC1$logFC>2,]

dim(diffB1) # 253   6
dim(diffB2)
dim(diffB3)
dim(diffC1)#145   6
write.csv(diffB1,'FOS_sig.csv')
write.csv(diffB2,'HSF_sig.csv')
write.csv(diffB3,'TP53_sig.csv')
write.csv(diffB4,'MYCN_sig.csv')
write.csv(diffC1,'MYC_ME_sig.csv')
