
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
expFile="rna-seq.txt"               
clusterFile="bulkRNA-group.txt"     
gmtFile="melanosome.gmt"                    

#1
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#2
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
param <- gsvaParam(exprData = data, geneSets = geneSets, minSize = 10, maxSize = 500)
gsva_result <- gsva(param, verbose = TRUE)
gsvaOut=rbind(id=colnames(gsva_result), gsva_result)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#3
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F,row.names = 1)

gsvaResult=gsva_result
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)


#####Single-cell GSVA analysis####
library(Seurat)
library(ggplot2)
library(clustree)
library(tidyverse)
library(cowplot)
library(dplyr)
library(patchwork)
library(UCell)
library(stringr)
library(ggplot2)
library(viridis)

####Scoring for datasets####
library(clusterProfiler)
library(GSVA)
getwd()
metabolism <- read.gmt("melanosome.gmt") 
unique(metabolism$term)
#Just repeat all the pathway names in the input file and run them one by one.
Oxidative <- subset(metabolism, term=="MELANOSOME TRANSPORT")
Oxidative <- list(Oxidative$gene)
#将基因整成list
names(Oxidative)[1] <- 'Oxidative'
load('melanocyte.RData')
DefaultAssay(sce) <- 'RNA'
metabolism_score <- AddModuleScore_UCell(sce,
                                         features=Oxidative, 
                                         name="_metabolism_score")
#耗时较久
#可视化所有细胞
pdf("melanogenesis-scRNA.pdf",height=4,width=4)
FeaturePlot(metabolism_score,features = "Oxidative_metabolism_score",            
            order = T,cols = viridis(256))
dev.off()
pdf("MELANOSOME TRANSPORT-scRNA.pdf",height=4,width=8)
FeaturePlot(metabolism_score,features = "Oxidative_metabolism_score",            
            order = T,cols = viridis(256), split.by = 'Group')
dev.off()
##箱线图可视化使用
library(ggrastr)
library(dplyr)
data<- FetchData(metabolism_score,vars = c("Group","Oxidative_metabolism_score"))
data$cellid <- case_when(data$seurat_clusters ==unique(data$seurat_clusters)[1] ~ "0",
                         data$seurat_clusters ==unique(data$seurat_clusters)[2] ~ '1',
                         data$seurat_clusters ==unique(data$seurat_clusters)[3] ~ '2',
                         data$seurat_clusters ==unique(data$seurat_clusters)[4] ~ '3', 
                         data$seurat_clusters ==unique(data$seurat_clusters)[5] ~ '4',
                         data$seurat_clusters ==unique(data$seurat_clusters)[6] ~ '5')

colors <- c('#507BA8','#F38D37','#5D9F53','#B5972D','#48998E','#E05758','#F1CE60')
colors <- c('#507BA8','#F38D37','#5D9F53','#B5972D','#48998E','#E05758')
colors <- c('#507BA8','#E05758')
pdf("MELANOSOME TRANSPORT-scRNA-boxplot.pdf",height=6,width=6)
ggplot(data, aes(x=Group,y=Oxidative_metabolism_score,fill=Group,color=Group)) +  
  theme_bw()+RotatedAxis()+  theme(panel.grid = element_blank(),        
                                   axis.text.x=element_text(size=12),        
                                   axis.text.y = element_text(size=10),       
                                   plot.title = element_text(hjust = 0.5),        
                                   legend.position = 'none')+  
  labs(x=NULL,y=NULL,title = "Oxidative_metabolism_score")+   
  geom_boxplot(position=position_dodge(0))+  scale_fill_manual(values = colors)+  
  geom_boxplot(position=position_dodge(0),color='black',               
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  comparisons =my_comparisons,
  label="p.signif",
  bracket.size=0.8, 
  size=6)
dev.off()

library(ggplot2)
library(ggpubr)
gene_score=metabolism_score
df<- FetchData(gene_score,vars = c("Group","gene_score"))
df<- FetchData(metabolism_score,vars = c("Group","Oxidative_metabolism_score"))
df$Group<- factor(df$Group,levels = c("NS","GCMN"))#设置顺序

#设置比较组
my_comparisons <- list(c("NS","GCMN")) 
#做小提琴图
pdf("GSVAscore_scRNA.pdf",height=6,width=6)
ggplot(df,aes(x=Group,y=Oxidative_metabolism_score,fill=Group))+  
  geom_violin(color='black',size=1)+theme_classic() +
  theme(text = element_text(size=10, colour = "black")) +   
  theme(plot.title = element_text(hjust = 0.5, size = 15),        
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 10),
  axis.title.y = element_text(color = 'black', size = 12),
  axis.line = element_line(size = 1))+   theme(legend.position="none") +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+stat_compare_means(method="t.test",hide.ns = F,
  comparisons =my_comparisons,
label="p.signif",
bracket.size=0.8, 
size=6)#添加显著性比较
dev.off() 

####TBX2 expression and GSVA score#####
mymatrix <- as.data.frame(metabolism_score@assays$RNA@data)
mymatrix<-t(mymatrix)%>%as.data.frame()
TBX2=as.data.frame(mymatrix$TBX2)

df1=as.data.frame(metabolism_score@meta.data$Oxidative_metabolism_score)
merged_df <- cbind(df1, TBX2)
colnames(merged_df) <- c("TBX2", "GSVAscore")

library(ggpubr)
p1 <- ggscatter(
  merged_df, 
  x = "TBX2", 
  y = "GSVAscore", 
  fill = "#4472C4",          # 散点填充色（匹配第一张图的蓝色）
  color = "#4472C4",         # 散点边框色（与填充色一致，使点更均匀）
  add = "reg.line",          # 添加回归线
  conf.int = TRUE,           # 显示置信区间
  add.params = list(
    color = "black",         # 回归线颜色（黑色，与第一张图一致）
    fill = "#B4C7E7"         # 置信区间填充色（浅蓝色，与散点协调）
  ),
  cor.coef = TRUE,           # 显示相关系数
  cor.coeff.args = list(
    method = "spearman",     # 改为 Spearman 相关（匹配第一张图）
    label.x = 0.1,           # 调整相关系数标签的 X 位置
    label.y = max(merged_df$GSVAscore, na.rm = TRUE) * 0.9  # 调整 Y 位置（避免与点重叠）
  ),
  cor.method = "spearman",   # 明确相关方法为 Spearman
  xlab = "Melanocyte proliferation score",  # X 轴标题（根据数据含义调整，如第一张图的 A/B 子图）
  ylab = "TBX2 Expression Level"             # Y 轴标题（匹配第一张图）
)

pdf(file="cor-TBX2 and GSVAscore.pdf", width=4, height=4)
print(p1)
dev.off()

#compare GSVAscore in high- and low-TBX2 group
tbX2_expr <- metabolism_score@assays$RNA@data["TBX2", ]  # 从稀疏矩阵提取一行（类型为 dgCMatrix 行）
tbX2_expr <- as.numeric(tbX2_expr)  # 转换为数值向量（关键修正！）
tbX2_mean <- mean(tbX2_expr, na.rm = TRUE)
high <- data[tbX2_expr > tbX2_mean, ]
low <- data[tbX2_expr <= tbX2_mean, ]

# 定义需要比较的评分列名
scores <- c(
  "melanocyte proliferation",
  "cell cycle and cell proliferation")

# 进行比较
results <- list()
for (score in scores) {
  res <- wilcox.test(high[[score]], low[[score]], alternative = "two.sided")
  results[[score]] <- res
}

# 查看结果
for (score in names(results)) {
  cat("Score:", score, "\n")
  print(results[[score]])
  cat("\n")
}

pathway_scores <- metabolism_score@meta.data[, grep("_score$", colnames(metabolism_score@meta.data))]
#merge
tbx2_pathway <- data.frame(
  TBX2 = tbx2_expr,
  pathway_scores
)

# 检查数据类型（确保均为数值型）
str(tbx2_pathway)

for (col in pathway_cols) {
  high_score <- na.omit(high_tbx2[[col]])
  low_score <- na.omit(low_tbx2[[col]])
  
  if (length(high_score) < 10 || length(low_score) < 10) {
    warning(sprintf("skip", col))
    next
  }
  
  plot_data <- data.frame(
    group = c(rep("High TBX2", length(high_score)), rep("Low TBX2", length(low_score))),
    score = c(high_score, low_score)
  )
  
  ggplot(plot_data, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = col, x = "TBX2 Group", y = "Pathway Score") +
    theme_bw()
  
  ggsave(sprintf("%s_TBX2_boxplot.pdf", col), width = 4, height = 4)
}
