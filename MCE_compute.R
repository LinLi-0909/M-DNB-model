###fpkm
fpkm <- read.csv('fpkm.csv',header = T,row.names = 1)
head(fpkm)
a <- rowMeans(log2(fpkm+1))
f_data <- fpkm[a>2,]
dim(f_data)
write.csv(f_data,'data_filter.csv')
entrez_id <- read.csv('entrez_id.csv',header = T)
head(entrez_id)
dim(entrez_id)
load('hprdAsigH-13Jun12.Rd')
entrez_id <- entrez_id[entrez_id$entrez_id%in%rownames(hprdAsigH.m),]
dim(entrez_id)
head(entrez_id)
f_data$gene <- rownames(f_data)
head(f_data)
data <- merge(entrez_id,f_data,all = F,sort = F)
dim(data) #5320   16
head(data[,1:5])
data1 <- data[,-1]
data2 <- aggregate(data1[,-1],by = list(data1$entrez_id),sum)
dim(data2) #5319   15
head(data2[,1:5])
colnames(data2)[1] <- 'entrez_id'
net <- hprdAsigH.m[rownames(hprdAsigH.m)%in%data2$entrez_id,rownames(hprdAsigH.m)%in%data2$entrez_id]
dim(net)
net_g <- data.frame(entrez_id=rownames(net))
data3 <- merge(net_g,data2,all = F,sort = F)
head(data3[,1:5])
rownames(data3) <- data3$entrez_id
head(data3[,1:5])
data3 <- data3[,-1]
data4 <- log2(data3+1.1)
head(data4[,1:5])
head(net[,1:5])
n <- ncol(net)
n
net1 <- net + diag(1,nrow = n,ncol = n)
dim(net1)
write.csv(net1,'MCE_net.csv')
write.csv(data4,'MCE_data.csv')

MCE <- readxl::read_xlsx('MCE.xlsx')
head(MCE)
MCE <- as.data.frame(MCE)
library(ggplot2)
head(MCE)
MCE$group <- factor(MCE$group,levels = c('ES','FOS','HSF1','MYCN','ME','MYC_ME','DE'))
ggplot(MCE, aes(x=group, y=MCE, color=group)) +scale_color_brewer(palette="Paired")+ 
  geom_line() + geom_point(size=10)+theme_bw()+ylab('MCE')
