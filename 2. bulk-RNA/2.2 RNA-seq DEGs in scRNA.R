folder_path <- "/home/Tanwei/GCMN/GCMN-7/bulkRNA"
# 使用dir.create创建文件夹
dir.create(folder_path)

experiment.merged <- readRDS("GCMN.rds")

####Draw a heatmap of the expression levels of transcriptome differential genes in single cells####
library(tibble)
library(tidyverse)
library(pheatmap)

markergene=read.table("type_marker_genes_GCMN.txt",
                      stringsAsFactors = F,
                      sep = "\t",
                      header = T,
                      fill = T,
                      comment.char = "#",
                      quote = "")

#downgenes.txt and upgene.txt comes from bulk-RNA data 'diffEXp.txt'
#logFC < '-2'-----downgene.txt
#logFC >2-----upgene.txt
rt1=read.table("downgenes.txt",
              stringsAsFactors = F,
              sep = "\t",
              header = T,
              fill = T,
              comment.char = "#",
              quote = "")
row.names(rt1)=rt1[,1]
same=intersect(rownames(rt1),as.vector(markergene[,7]))
row.names(markergene)=markergene[,7]
markergene=as.matrix(markergene)
markergene1=markergene[same,]

rt2=read.table("upgenes.txt",
               stringsAsFactors = F,
               sep = "\t",
               header = T,
               fill = T,
               comment.char = "#",
               quote = "")
row.names(rt2)=rt2[,1]
same=intersect(rownames(rt2),as.vector(markergene[,7]))
row.names(markergene)=markergene[,7]
markergene=as.matrix(markergene)
markergene2=markergene[same,]

#将数据储存为RData
save(experiment.merged,markergene1,file="GCMN-heatmap-down.RData")
save(experiment.merged,markergene2,file="GCMN-heatmap-up.RData")
load("GCMN-heatmap-down.RData")
library(scRNAtoolVis)
p2 <- averageHeatmap(object = experiment.merged1,
                     markerGene = markergene1, 
                     group.by = "celltype")
pdf("rnaseq-heatmap-down.pdf", width = 12, height = 12)
# 绘制热图
draw(p2 )
# 关闭图形设备
dev.off()

library(scRNAtoolVis)
p2 <- averageHeatmap(object = experiment.merged1,
                     markerGene = markergene2, 
                     group.by = "celltype")
pdf("rnaseq-heatmap-up.pdf", width = 12, height = 12)
# 绘制热图
draw(p2 )
# 关闭图形设备
dev.off()
