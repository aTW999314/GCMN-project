###细胞群的相关性分析###
library(Seurat)

install.packages("corrplot")
install.packages("ggplot2")
install.packages("ggpubr")
library(corrplot)
library(ggplot2)
library(ggpubr)

###两两相关性点图绘制####
#学习资源见https://zhuanlan.zhihu.com/p/586234818
TBX2=as.data.frame(mymtrix$TBX2)

cor.test (gene3, gene9, method="spearman")
#显示结果
#Spearman's rank correlation rho

#data:  gene3 and gene9
#S = 198320, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.6677523 

library(ggplot2)
#提取单细胞矩阵中某个基因的表达
install.packages("data.table")
library(data.table)
#利用system.time记录运行时间
system.time({fwrite(x = as.data.frame(uterus[["RNA"]]@counts), row.names=T,file = "counts.csv")})
rt <- read.csv("counts.csv", header = TRUE)
df<- rt[, "TBX2"]


###提取单细胞数据的表达矩阵###
mymatrix <- as.data.frame(sce@assays$RNA@data)
mymatrix<-t(mymatrix)%>%as.data.frame()
TBX2=as.data.frame(mymatrix$TBX2)

# 假设df1和df2的行数相同
df1=as.data.frame(metabolism_score@meta.data$Oxidative_metabolism_score)
merged_df <- cbind(df1, TBX2)
colnames(merged_df) <- c("TBX2", "senescence")



# 提取名为"Gene2"的基因表达矩阵
gene_expression <- subset(df, Gene == "Gene2")

p1 <- ggplot(merged_df, aes(x = metabolism_score@meta.data$Oxidative_metabolism_score, y = mymatrix$TBX2))
p2 <- p1 + geom_point() 
p3 <- p2 + geom_smooth(method="lm")
p3

#相关性图形
length=length(levels(factor(data$geneCluster)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(m6Ascore, TMB)) + 
  xlab("m6Ascore")+ylab("Tumor Burden Mutation")+
  geom_point(aes(colour=geneCluster))+
  scale_color_manual(values=bioCol[1:length])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =m6Ascore, y =TMB))
#相关性图形
pdf(file="cor-senescence.pdf", width=4, height=4)
print(p3)
dev.off()

####添加相关性系数与p值####
install.packages("ggpubr")
library(ggpubr)
p1=ggscatter(merged_df, x = "TBX2", y = "senescence", 
          fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson")

pdf(file="cor-pyroptosis.pdf", width=4, height=4)
print(p1)
dev.off()






