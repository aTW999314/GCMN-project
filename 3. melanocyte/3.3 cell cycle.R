###cell cycle analysis###
library(Seurat)
library(tidyverse)
library(ggsci)
load('melanocyte.RData')
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
test.seu <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

head(test.seu@meta.data,2)
pdf("cellcycle1.pdf",width = 10,height = 4)
test.seu@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
dev.off()

pdf("cellcycle2.pdf",width = 5,height = 4)
DimPlot(test.seu,reduction = "umap")
dev.off()

#提取细胞行名以及cluster分类以及细胞的周期评分信息
cell_barcodes <- rownames(test.seu@meta.data)  # 直接获取行名作为细胞名
cluster_info <- test.seu@meta.data[, "seurat_clusters"]  # 假设列名为 "seurat_clusters"
cell_cycle_scores <- test.seu@meta.data[, c("S.Score", "G2M.Score")] %>%
  dplyr::mutate(G1.Score = 1 - S.Score - G2M.Score)  # 计算 G1 期评分
library(dplyr)
# 将行名（细胞名）转换为列
combined_data <- data.frame(
  cell_barcode = cell_barcodes,
  cluster = cluster_info,
  cell_cycle_scores
)

####保存数据####
save(combined_data,test.seu,file="seu.Rdata")
head(combined_data)
#1. 整体差异检验（多组比较）
#若有 ≥3 个 Cluster，先通过 Kruskal-Wallis 检验 判断所有组间是否存在显著差异（适用于非正态分布数据）：
# 检验假设：所有Cluster的S期评分中位数相等
kruskal_result <- kruskal.test(G1.Score ~ cluster, data = combined_data)
print(kruskal_result)
# 输出结果包含H统计量和p值（若p < 0.05，说明至少有两组存在差异）

###2. 两两比较（多重检验校正）
#若整体检验显著，进一步进行 两两 Wilcoxon 秩和检验，并通过 Benjamini-Hochberg 法 校正 FDR（控制假发现率）：
library(dplyr)
library(purrr)
# 获取所有Cluster对
cluster_pairs <- combn(unique(combined_data$cluster), 2, simplify = FALSE)

# 定义两两检验函数
pairwise_wilcox <- function(pair) {
  group1 <- combined_data$G1.Score[combined_data$cluster == pair[1]]
  group2 <- combined_data$G1.Score[combined_data$cluster == pair[2]]
  wilcox_result <- wilcox.test(group1, group2, alternative = "two.sided")
  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = wilcox_result$p.value
  )
}

# 执行所有两两检验并整理结果
pairwise_results <- map_dfr(cluster_pairs, pairwise_wilcox)

# 计算FDR校正后的p值（q值）
pairwise_results <- pairwise_results %>%
  mutate(q_value = p.adjust(p_value, method = "fdr"))

print(pairwise_results)
# 输出每对Cluster的p值和q值，q < 0.05视为显著差异


library(ggplot2)
pdf("G1-cluster score.pdf",width = 6,height = 6)
ggplot(combined_data, aes(x = cluster, y = G1.Score, color = cluster)) +
  # 箱线图：展示中位数、四分位数范围
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +  
  # 点图：每个点代表一个细胞，添加抖动避免重叠
  geom_jitter(shape = 21, fill = "white", size = 2, width = 0.2, alpha = 0.8) +  
  # 美化：删除颜色图例（与x轴标签重复），添加标题和轴标签
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Cluster", y = "G1.Score", title = "G1 score")

dev.off()

install.packages("ggpubr")
library(ggpubr)
pdf("S-cluster score with pvalue.pdf",width = 6,height = 6)
ggplot(combined_data, aes(x = cluster, y = S.Score, color = cluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) +  
  geom_jitter(shape = 21, fill = "white", size = 2, width = 0.2, alpha = 0.8) +  
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Cluster", y = "S.Score") +
  (*: p < 0.05, **: p < 0.01, ***: p < 0.001）
   stat_compare_means(
     method = "wilcox.test",          # 使用Wilcoxon检验
     label = "p.signif",              # 显示星号
     comparisons = cluster_pairs       # 指定比较的Cluster对
   )
dev.off()
   
# 计算每个 cluster 的 S 期评分平均值和方差
cluster_stats <- combined_data %>%
     group_by(cluster) %>%  # 按 cluster 分组
     summarise(
       mean_G2M = mean(G2M.Score, na.rm = TRUE),  # 平均值（na.rm=TRUE 忽略缺失值）
       var_G2M = var(G2M.Score, na.rm = TRUE)      # 方差（样本方差，分母为 n-1）
     )
   
 # 查看结果
print(cluster_stats)


all_genes <- rownames(seu@assays$RNA@data)

# 检查目标基因是否存在（返回 TRUE/FALSE）
target_genes <- c("TBX2")
target_genes %in% all_genes

# 提取数据（假设 cluster 列名为 "seurat_clusters"）
expression_data <- FetchData(
  test.seu,
  vars = c(target_genes, "seurat_clusters","Phase","Group","orig.ident","S.Score","G2M.Score")  # vars 参数指定要提取的基因和 meta.data 列
)

write_xlsx(expression_data, path = "TBX2 and gene expression of cycle and malanogenesis.xlsx")

view(g2m.genes)
view(s.genes)

write_xlsx(kruskal_results, path = "kruskal_results.xlsx")

#TBX2 expression in phase
expression_data$Phase <- factor(expression_data$Phase, levels = c("G1", "S", "G2M"))

pdf("TBX2_in_phase.pdf", width = 8, height = 6)
p <- ggplot(expression_data, aes(x = Phase, y = TBX2, fill = Phase)) +
  geom_boxplot(alpha = 0.8) +  # 箱图
  geom_jitter(                  # 抖动点（展示原始数据分布）
    shape = 21, fill = "white", size = 1, 
    width = 0.2, alpha = 0.4
  ) +
  theme_minimal() +            # 简洁主题
  labs(                        # 图表标签
    x = "Cell Cycle Phase", 
    y = "TBX2 Expression", 
    title = "TBX2 Expression Across Cell Cycle Phases"
  ) +
  theme(                       # 调整文本格式
    legend.position = "none",  # 隐藏图例
    axis.text.x = element_text(angle = 45, hjust = 1),  # x轴标签倾斜
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加Kruskal-Wallis检验结果（多组差异检验）
  stat_compare_means(
    method = "kruskal.test",  # 非参数多组检验（适用于非正态分布数据）
    label = "p.format",       # 显示p值（如"p = 0.001"）
    label.y = max(expression_data$TBX2, na.rm = TRUE) * 1.1  # p值位置（y轴顶部）
  )

print(p)
dev.off()

# 统计各cluster的周期阶段比例
cluster_phase <- test.seu@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(phase_proportion = cell_count / sum(cell_count)) %>%
  ungroup()


# 可视化：各cluster的周期阶段比例（堆叠条形图）
pdf("cluster in phase.pdf",width = 6,height = 6)
ggplot(cluster_phase, aes(x = factor(seurat_clusters), y = phase_proportion, fill = Phase)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster", y = "周期阶段比例", fill = "周期阶段") +
  theme_minimal() +
  ggtitle("各Cluster的细胞周期阶段分布")
dev.off()

