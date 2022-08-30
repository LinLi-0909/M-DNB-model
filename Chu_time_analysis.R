library(Seurat)
time <- readRDS('time.rds')
time
#15914 features across 758 samples within 1 assay 
time <- NormalizeData(time, normalization.method = "LogNormalize", scale.factor = 1000000)
all.genes <- rownames(time)
time <- FindVariableFeatures(time,selection.method = "vst",nfeatures = 5000)
time <- RunPCA(time, features =VariableFeatures(time))
ElbowPlot(time)
time <- RunTSNE(time, dims = 1:20)
DimPlot(time, reduction = "tsne",group.by = 'cluster',label.size = 4,pt.size = 2)

