library(future) # parrel for Seurat
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(viridis)
library(scales)
library(scRNAtoolVis)

####Extract melanocyte for further in-depth analysis####
View(experiment.merged@meta.data$celltype)
table(experiment.merged@meta.data$Group)
sce = experiment.merged[,experiment.merged@meta.data$celltype %in% c('melanocytes')]
#also can be Keratinocyte cells or other cells

#FindVariableFeatures
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst",
                            nfeatures = 2000)

#Show the overall expression level after standardization
top10 <- head(x = VariableFeatures(sce), 10)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf("Melanocyte-Norm-feature_variable_plot.pdf",width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
#### Normalization and PCA####
#Normalization
sce <- ScaleData(
  object = sce,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.mt"))

#PCA
sce <- RunPCA(object = sce, 
              features = VariableFeatures(Melanocyte),
              verbose = F,npcs = 20)

#PCA Results Display-1
pdf("Melanocyte-PCA-VizDimLoadings.pdf",width = 7,height = 5)
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
dev.off()

#PCA Results Display-2
pdf("Melanocyte-PCA-DimPlot.pdf",width = 5,height = 4)
DimPlot(sce, reduction = "pca")
dev.off()

#PCA Results Display-3
pdf("Melanocyte-PCA-DimHeatmap.pdf",width = 5,height = 4)
DimHeatmap(sce, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

sce <- JackStraw(sce, num.replicate = 100,dims = 40)
sce <- ScoreJackStraw(sce, dims = 1:20)
pdf("Melanocyte-PCA-JackStrawPlot_40.pdf",width = 6,height = 5)
JackStrawPlot(object = sce, dims = 1:20)
dev.off()

pdf("Melanocyte-PCA-ElbowPlot.pdf",width = 6,height = 5)
ElbowPlot(sce,ndims = 20)
dev.off()

dim.use <- 1:12
#### Cell clustering TSNE algorithm ####
#TSNE algorithm
sce <- FindNeighbors(sce, dims = dim.use)
sce <- FindClusters(sce, resolution = 0.1)
table(sce@meta.data$celltype)

sce <- RunTSNE(sce, dims = dim.use, 
               do.fast = TRUE)
pdf("Melanocyte-TSNEPlot_res0.1_PC.pdf",width = 5,height = 4)
DimPlot(object = sce, pt.size=0.5,label = T)
dev.off()

#uAMP
sce <- RunUMAP(sce, dims = dim.use, 
               do.fast = TRUE)
pdf("Melanocyte-CellCluster-UMAPPlot_res0.1-PC.pdf",width = 5,height = 4)
DimPlot(object = sce, pt.size=0.5,label = T, reduction = "umap")
dev.off()

Group <- 'GCMN'
sce1 <- subset(sce, Group=='GCMN')
sce2 <- subset(sce, Group=='NS')
sce3=merge(sce1,sce2)

pdf("Melanocyte-CellCluster-umaplot_group-PC1.pdf",width = 5,height = 4)
DimPlot(object = sce3, 
        group.by="Group", 
        pt.size=0.5,reduction = "umap")
dev.off()

#cell annotation
sce.markers1 <- FindAllMarkers(sce, only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.3)

write.table(sce.markers1,
            file="GCMN_Melanocyte_marker_gene_res0.1.txt",
            sep="\t",quote = F,row.names = F)


#Draw a cell proportion figure
#Calculate the proportion of different cells in each subcluster
meta_data <- sce@meta.data 
plot_data <- data.frame(table(meta_data$Group,meta_data$seurat_clusters))
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
pdf("Melanocytes-propotion-group.pdf",width = 8,height = 4)
ggplot(plot_data,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + labs(fill = "Cluster")
dev.off()

Cellratio <- prop.table(table(sce$Group,Idents(sce)), margin = 2)
#Calculate the proportion of different cell populations in each group of samples
Cellratio <- as.data.frame(Cellratio)
table(sce$Group)


#top-marker gene dotplot
rt <- c("FABP7","TFAP2B","MCOLN3","SPP1","KCNJ13","TBX2",
        "VMP1","SQSTM1","PPP1R10","SERPINB4","AQP3","SFN","PHACTR1",
        "TNFRSF12A","GADD45B","CRABP1","PCSK2","SLC6A17")

sce$celltype <- Idents(sce)
head(sce@meta.data)
pdf("MarkerGene-melanocyte_Dotplot.pdf",width = 13,height = 5)
DotPlot(sce, features = rt)+
  RotatedAxis()
dev.off()