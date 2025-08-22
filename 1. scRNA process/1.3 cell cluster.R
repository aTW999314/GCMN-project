library(future) # parrel for Seurat
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(viridis)
library(scales)
library(scRNAtoolVis)

dim.use <- 1:15
####Cell clustering TSNE algorithm ####

sce <- FindNeighbors(sce, dims = dim.use)
sce <- FindClusters(sce, resolution = 0.1)

#RunTSNE
sce <- RunTSNE(sce, dims = dim.use, 
               do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.1_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, pt.size=0.5,label = T)
dev.off()

#RunuAMP
sce <- RunUMAP(sce, dims = dim.use, 
               do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-UMAPPlot_res0.1_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, pt.size=0.5,label = T, reduction = "umap")
dev.off()

#orig.ident
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 5,height = 4)
pdf("gcmn-umap2.pdf",width = 8,height = 4)
DimPlot(object = sce, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "umap")
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 8,height = 4)
DimPlot(object = sce, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

#### Calculate marker genes ####
all.markers <- FindAllMarkers(sce, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 2)

write.table(all.markers,
            file=paste0("GCMN_cluster_marker_gene_","PC.txt"),
            sep="\t",quote = F,row.names = F)

#####11. Annotate single-cell data####
#Manually annotate by reviewing literature based on differentially expressed genes
sce <- RenameIdents(
  object = sce,
  "0" = "melanocytes","1" = "Keratinocyte",
  "2" = "Keratinocyte",
  "3" = "Fibroblast",
  "4" = "Endothelial cells", "5" = "Melanocytes","6" = "melanocytes",
  "7" = "Fibroblasts","8" = "melanocytes","9" = "T/DC",
  "10" = "Myeloid cells","11" = "Fibroblast","12" = "melanocytes",
  "13" = "Mast cells","14" = "Endothelial cells")

table(sce$RNA_snn_res.0.1)
pdf("cell_type_GCMN-7.pdf",width = 7,height = 4)
DimPlot(object = sce, pt.size=0.5,
        reduction = "umap",label = T) +
  ggsci::scale_color_igv()
dev.off()

df <- sce@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = sce@meta.data$celltype)
colnames(df)

library(ggplot2)
library(dittoSeq)
library(ggrastr)
library(tidydr)
library(scatterpie)

#Calculate the cell proportion
Cellratio1 <- prop.table(table(Idents(sce),sce$orig.ident), margin = 2)
#Calculate the proportion of different cell populations in each group of samples
Cellratio1 <- as.data.frame(Cellratio1)
table(sce$Group)
write.table( Cellratio1, file='Cellratio.txt',sep = "\t")
library(tidyr)
freq <-spread(Cellratio, Var1, Freq)
colnames(freq)[1] <- 'celltype'
freq <- freq[sort(freq$celltype),]

label <- df %>%group_by(cell_type) %>%
  summarise(umap_1 = median(umap_1),
            umap_2 = median(umap_2))%>%
  as.data.frame()
rownames(label) <- label$cell_type
label <- label[c(5,3,7,6,1,2,4), ]


cell_number <- as.data.frame(table(sce$celltype))
cell_number <- cell_number[c(5,3,7,6,1,2,4), ]
colnames(cell_number)[2]<-'cellnumber'
cell_number$cellnumber <- log2(cell_number$cellnumber)/3

data = cbind(freq,label[,c(2:3)], cell_number[,c(2)])
colnames(data)[6]<- 'cellnumber'

#提取GCMN和NS的先进行后续分析
Group <- 'GCMN'
sce1 <- subset(experiment.merged1, Group=='GCMN')
Group <- 'NS'
sce2 <- subset(experiment.merged1, Group=='NS')
sce=merge(sce1,sce2)
table(sce$Group)


pdf("cell_type_GCMN.pdf",width = 7,height = 4)
DimPlot(object = sce1, pt.size=0.5,
        reduction = "umap",label = T) +
  ggsci::scale_color_igv()
dev.off()

pdf("cell_type_NS.pdf",width = 7,height = 4)
DimPlot(object = sce2, pt.size=0.5,
        reduction = "umap",label = T) +
  ggsci::scale_color_igv()
dev.off()

#Add Celltype information in the metadata
sce$celltype <- Idents(sce)
sce$seurat_clusters=sce$celltype
sce
head(sce@meta.data)
str(all.markers)

#Draw a cell proportion chart
#Calculate the proportion of different cells in each type of sample and plot it.
meta_data <- sce@meta.data 
table(sce@meta.data$orig.ident)

plot_data1 <- data.frame(table(meta_data$Group,meta_data$seurat_clusters))
plot_data1$Total <- apply(plot_data1,1,function(x)sum(plot_data1[plot_data1$Var1 == x[1],3]))
plot_data1 <- plot_data1 %>% mutate(Percentage = round(Freq/Total,3) * 100)

pdf("group_melanocyte.pdf",width = 4,height = 4)
ggplot(plot_data1,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + labs(fill = "seurat_clusters")
dev.off()
getwd()

#Compare marker genes again based on cell annotation
all.markers1 <- FindAllMarkers(sce, only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.3)
write.table(all.markers1,
            file="type_marker_genes_GCMN.txt",
            sep="\t",quote = F,row.names = F)

#Draw a cell proportion chart
saveRDS(sce,"GCMN.rds") 
experiment.merged <- readRDS("GCMN.rds")
count_table <- table(experiment.merged@meta.data$celltype, experiment.merged@meta.data$orig.ident)
count_table
install.packages('ggstatsplot')
rt=as.data.frame(count_table)

write.table(rt,file="count_table.txt",sep="\t",row.names=F,quote=F)

pdf('Cell_prop1.pdf', width = 12, height = 7)
rt %>% 
  ggbarstats(x = Var2, 
             y = Var1,
             counts = Freq)+
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45))
dev.off()


#DoHeatmap
#can revise，features=c(“gene1”,"gene2","gene3","gene4","gene5")
rt=c("MLANA","DCT","CDH19","KRT5","KRT14","KRT17","COL1A1","COL1A2","COL6A2",
     "CD74","FLT1","PECAM1","PTPRC","SRGN","SAMSN1",
     "HLA-DRA","HLA-DRB1","LYZ","TPSB2","TPSAB1","CPA3")

#top-marker gene dotplot
pdf("MarkerGene-celltype_umap1.pdf",width = 13,height = 5)
DotPlot(sce, features = rt)+
  RotatedAxis()
dev.off()
