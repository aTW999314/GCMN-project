.libPaths()
install.packages("sp")
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(stringi)
library(stringr)

GCMN_S1 <- Read10X("./GCMN_S1/")
colnames(GCMN_S1) <- paste(colnames(GCMN_S1),"GCMN_S1",sep = "_")

GCMN_S2 <- Read10X("./GCMN_S2/")
colnames(GCMN_S2) <- paste(colnames(GCMN_S2),"GCMN_S2",sep = "_")

NS_S1<- Read10X("./NS_S1/")
colnames(NS_S3) <- paste(colnames(NS_S3),"NS_S3",sep = "_")

NS_S2<- Read10X("./NS_S2/")
colnames(NS_S4) <- paste(colnames(NS_S4),"NS_S4",sep = "_")

GCMN_S3<- Read10X("./S2/")
colnames(GCMN_S3) <- paste(colnames(GCMN_S5),"GCMN_S3",sep = "_")

NS_S3<- Read10X("./S3/")
colnames(NS_S6) <- paste(colnames(NS_S6),"NS_S6",sep = "_")

GCMN_S1.df <- data.frame(GCMN_S1) %>% mutate(Gene = rownames(GCMN_S1))
GCMN_S2.df <- data.frame(GCMN_S2) %>% mutate(Gene = rownames(GCMN_S2))
GCMN_S3.df <- data.frame(GCMN_S3) %>% mutate(Gene = rownames(GCMN_S2))

NS_S1.df <- data.frame(NS_S1) %>% mutate(Gene = rownames(NS_S1))
NS_S2.df <- data.frame(NS_S2) %>% mutate(Gene = rownames(NS_S2))
NS_S3.df <- data.frame(NS_S3) %>% mutate(Gene = rownames(NS_S3))

experiment.data <- merge(GCMN_S1.df, GCMN_S2.df,by = "Gene", all.x = TRUE)
experiment.data <- merge(experiment.data, GCMN_S3.df,by = "Gene", all.x = TRUE)
experiment.data <- merge(experiment.data, NS_S1.df, by = "Gene", all.x = TRUE)
experiment.data <- merge(experiment.data, NS_S2.df, by = "Gene", all.x = TRUE)
experiment.data <- merge(experiment.data, NS_S3.df, by = "Gene", all.x = TRUE)

rownames(experiment.data) <- experiment.data$Gene
experiment.data <- experiment.data %>% select(-Gene)
experiment.data[is.na(experiment.data)] <- 0

sam.name <- "GCMN"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi", 
  min.cells = 20,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

str(experiment.aggregate)
save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_SeuratObject.RData"))

dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)
table(sce@meta.data[["orig.ident"]])
table(sce@meta.data[["orig.ident"]])

experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

pdf("QC_UMI.pdf",width = 5,height = 4.5)
sce@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = Group, fill=Group)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

pdf("QC_ncount_RNA.pdf",width = 5,height = 4.5)
hist(sce@meta.data$nFeature_RNA, breaks = 60)
abline(v = 2753, col = "red")
dev.off()

plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
pdf("ncount-RNA.pdf",width = 8,height = 4.5)
hist(experiment.aggregate@meta.data$nCount_RNA, main = "Histogram of Values", xlab = "nCount_RNA", freq = TRUE, ylab = "Frequency")
dev.off()
pdf("nfeature-RNA.pdf",width = 8,height = 4.5)
hist(experiment.aggregate@meta.data$nFeature_RNA, main = "Histogram of Values", xlab = "nCount_RNA", freq = TRUE, ylab = "Frequency")
dev.off()

experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nFeature_RNA < 6000 & 
                                 nCount_RNA > 1000 &
                                 percent.mt < 15)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")


experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 2000)

top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

#####filter2####
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.mt"))

experiment.aggregate <- RunPCA(object = sce, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 20)

pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 5,height = 4)
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

experiment.aggregate <- JackStraw(experiment.aggregate, num.replicate = 100,dims = 40)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:20)
pdf(paste0("./",sam.name,"/PCA-JackStrawPlot_40.pdf"),width = 6,height = 5)
JackStrawPlot(object = experiment.aggregate, dims = 1:20)
dev.off()

#碎石图
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(experiment.aggregate,ndims = 20)
dev.off()

dims <- 15
# 扫描 pK 取值，找到最优值（需耗时，可跳过直接用默认值）
sweep.res <- paramSweep_v3(
  experiment.aggregate,   # Seurat 对象
  PCs = 1:15,             # PCA 维度（与 RunPCA 一致）
  sct = FALSE             # 若用 SCTransform 归一化，设为 TRUE
)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal.pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]  # 最优 pK

# 设定预期双细胞数量（nExp）
nCells <- ncol(experiment.aggregate)  # 细胞总数
nExp <- round(nCells * 0.05)          # 假设 5% 双细胞比例（可调整）

# 运行 DoubletFinder
experiment.aggregate <- doubletFinder_v3(
  experiment.aggregate, 
  PCs = 1:10,           # PCA 维度
  pN = 0.25,            # 默认参数（无需修改）
  pK = optimal.pK,      # 若未优化，可用默认值（如 0.01）
  nExp = nExp,          # 预期双细胞数
  reuse.pANN = FALSE,   # 是否复用 pANN 矩阵
  sct = FALSE           # 归一化方法适配
)

# 统计双细胞比例
table(experiment.aggregate$DF.classifications)

# 过滤双细胞
experiment.aggregate <- experiment.aggregate[, 
                                             experiment.aggregate$DF.classifications == "Singlet"
]

#integration by harmony
sce=experiment.aggregate
pdf("umap_seurat_before_clusters.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="seurat_clusters")
dev.off()
pdf("umap_before_groups.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="Group")
dev.off()
pdf("umap_before_samples.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="orig.ident")
dev.off()
pdf("umap_before_batch.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="Batch")
dev.off()


##去批次数据分析
library(harmony)
#可以根据Batch或者orig.ident去批次
sce <- RunHarmony(sce, split.by , plot_convergence = F,dims.use = 1:dims)
sce = FindNeighbors(sce, reduction = "harmony",dims = 1:dims)
sce = FindClusters(sce, resolution = 0.1)
sce= RunTSNE(sce, reduction = "harmony", dims = 1:dims)
sce = RunUMAP(sce, reduction = "harmony", dims = 1:dims)
#去批次后绘图
pdf("after_seurat_clusters.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="seurat_clusters")
dev.off()
pdf("umap_after_groups.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="Group")
dev.off()
pdf("umap_after_samples.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="orig.ident")
dev.off()
pdf("umap_after_batch.pdf",height=7,width=10)
DimPlot(sce,reduction = "umap",group.by="Batch")
dev.off()

dim(sce)
