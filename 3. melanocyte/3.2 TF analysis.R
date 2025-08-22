#####TF analysis####
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)
library(ComplexHeatmap)

BiocManager::install(c("AUCell", "RcisTarget"),version='3.18')
BiocManager::install(c("GENIE3"))
BiocManager::install(c("zoo", "mixtools", "rbokeh"))

##1, For human:
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")
#can download directly using Linux, which is faster and more likely to succeed than R.
#wget  https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather
# mc9nr: Motif collection version 9: 24k motifs

dir.create("cisTarget_databases");   #Create a folder to save the database
setwd("cisTarget_databases")
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)

##Set up the analysis environment
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

mydbDIR <- "/home/Tanwei/GCMN/GCMN-7/Melanocyte/TF"
names(mydbs) <- c("500bp", "10kb")
mydbs <- c("hg19__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=5,
                                  dbDir=mydbDIR,dbs = mydbs)                                  )

saveRDS(scenicOptions, "int/scenicOptions.rds")

####Seurat => loom####
workPath <- '/home/Tanwei/GCMN/GCMN-7/Melanocyte/TF'
setwd(workPath)

####Prepare cell meta information,cellInfo####
scRNA=sce
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")

##Prepare the expression matrix
#To save computing resources, randomly select a subset of data from 1000 cells
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
scRNAsub <- scRNA
saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(sce@assays$RNA@counts)

####Shared Expression Network Computing####
##==Inference of transcriptional regulatory networks==##
##Gene filtering
#The filtering criteria are that the sum of gene expression levels is greater than 3% of the number of cells, and the gene is expressed in 1% of the cells.
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##Calculate the correlation matrix
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets correlation regression analysis
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#This step consumes a great deal of computing resources, and a personal computer would require several hours of running time.
save(exprMat_filtered_log, scenicOptions,file="genie3.RData")
load("genie3.RData")

##推断共表达模块
scenicOptions <-runSCENIC_1_coexNetwork2modules(scenicOptions)
save(scenicOptions,file="runSCENIC_1.RData")
#主要运行结果是int目录下的1.6_tfModules_asDF.Rds
#Target    TF method corr
#1   BCDIN3D ABCF2   w001    1
#2    STRADB ABCF2   w001    1
#3 HIST1H2BG ABCF2   w001    1
#4      GID4 ABCF2   w001    1
#5    LPCAT1 ABCF2   w001   -1
#6      ZHX2 ABCF2   w001    1
    
scenicOptions <-runSCENIC_2_createRegulons(scenicOptions)
save(scenicOptions,file="runSCENIC_2.RData")

scenicOptions <-runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
save(scenicOptions,exprMat,file="runSCENIC_3.RData")
load("runSCENIC_3.RData")

scenicOptions@fileNames$output["loomFile",] <- "output/GCMN.loom"
export2loom(scenicOptions, exprMat)

scenicOptions <-runSCENIC_4_aucell_binarize(scenicOptions)
save(scenicOptions,exprMat_filtered_log, file = "runSCENIC_4_aucell_binarize.RData")

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="TBX2"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="TBX2" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

####Visualization of SCENIC results using Seurat#### 
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')
scRNAauc=readRDS("scRNAauc.rds")

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

dir.create('scenic_seurat')
#FeaturePlot
p1 = FeaturePlot(scRNAauc, features='TBX2_extended_13g', label=T, reduction = 'umap')
ggsave('TBX2.pdf', p1, width=4 ,height=5)

#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "TBX2_extended_13g", group.by="seurat_clusters") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features = "TBX2_extended_13g", pt.size = 0, group.by="seurat_clusters") + 
  theme(legend.position='none')
plotc = p2
ggsave('scenic_seurat/Ridge-Vln_TBX2_extended_13g.pdf', plotc, width=8, height=8)

####regulon Expression heatmap####
library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'cluster')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#Select some interesting regulons
my.regulons <- c('PBX1_33g','PBX1_extended_36g','TBX2_extended_13g','SOX9_extended_11g','SIRT6_extended_43g',
                 'RORA_15g','THRA_extended_13g','ZNF471_extended_170g','NR2F2_extended_12g',
                 'HIVEP1_extended_14g','IKZF2_15g','EZH2_20g',
                 'NFKB2_extended_44g','RARA_extended_30g','NFKB2_37g','PBX3_25g','ZBTB7B_extended_93g',
                 'CLOCK_10g','PATZ1_extended_57g','MAFK_17g','BCL3_extended_26g')
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]

subcell <- sample(colnames(myAUCmatrix),1000)
myAUCmatrix1 <- myAUCmatrix[,subcell]
pdf(file = "cluster-myAUCmatrix_heatmap-RSS.pdf",width = 8,height = 8)
hm=pheatmap(myAUCmatrix1, show_colnames=F, annotation_col=celltype)
dev.off()

pdf(file = "cluster-myBINmatrix_heatmap-RSS.pdf",width = 8,height = 8)
hm=pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
            color = colorRampPalette(colors = c("white","black"))(100))
dev.off()

####celltype TF####
library(SCopeLoomR)
# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/GCMN.loom"
export2loom(scenicOptions, exprMat)
loom <- open_loom("./output/GCMN.loom")
# Read information from loom file:
regulons_incidMat<-get_regulons(loom,column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
regulons<-regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
embeddings

sub_regulonAUC<-regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)
identical(colnames(sub_regulonAUC),colnames(sce))

cellTypes<-data.frame(row.names=colnames(sce),
                      celltype=sce$celltype)
head(cellTypes)
sub_regulonAUC[1:4,1:4]
save(sub_regulonAUC,
     cellTypes,
     sce,
     file='for_rss_and_visual.Rdata')
table(sce$celltype)

selectedResolution<-"celltype"
cellsPerGroup<-split(rownames(cellTypes),cellTypes[,selectedResolution])

sub_regulonAUC<-sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)

#Regulon Specificity Score, RSS
selectedResolution<-"celltype"
rss<-calcRSS(AUC=getAUC(sub_regulonAUC),
             cellAnnotation=cellTypes[colnames(sub_regulonAUC),selectedResolution])
rss=na.omit(rss)
rssPlot<-plotRSS(rss,
                 labelsToDiscard=NULL,
                 zThreshold=1,
                 cluster_columns=FALSE,
                 order_rows=T,
                 thr=0.01,
                 varName="cellType",
                 col.low='#330066',
                 col.mid='#66CC66',
                 col.high='#FFCC33',
                 revCol=F,
                 verbose=TRUE
)

pdf(file = "cluster-specific-TF-heatmap-RSS1.pdf",width = 8,height = 10)
print(rssPlot1$plot)
dev.off()
rt=rssPlot$plot
write.table(rt,file="RSS_cluster.txt",sep="\t",row.names = TRUE)
#2.rank figure
plotRSS_oneSet(rss, setName = "0")
plotRSS_oneSet(rss, setName = "1")
plotRSS_oneSet(rss, setName = "5")

regulonActivity_byGroup<-sapply(cellsPerGroup,
                                function(x)
                                  rowMeans(getAUC(sub_regulonAUC)[,x]))
range(regulonActivity_byGroup)

regulonActivity_byGroup_Scaled<-t(scale(t(regulonActivity_byGroup),
                                        center=T,scale=T))
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

save(regulonActivity_byCellType_Scaled,file = "regulonActivity_byCellType_Scaled.RData")
regulons <- loadInt(scenicOptions, "regulons")

####Select the regulons with high rankings that you are interested in for heatmap display.
top5regulon = apply(regulonActivity_byCellType_Scaled, 2, function(x) names(x[order(abs(x))[1:5]]))
top5regulonV=top5regulon %>%as.vector()%>%unique()

pdf(file = "cluster-specific-TF-heatmap.pdf",width = 8,height = 14)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[top5regulonV,], name="Regulon activity", row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled[top5regulonV,])[row_order(hm)] # to save the clustered regulons for later
dev.off()

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$cluster),function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pdf(file = "cluster-specific-TF-heatmap-binary.pdf",width = 8,height = 14)
hm=pheatmap::pheatmap(binaryActPerc_subset,color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

#Heatmap of transcription factor expression levels
top3tfgene <- c("PBX1","TBX2","BHLHE41","ZNF471","NR2F2","IKZF2",
                "PBX3","RORA","SOX9","CLOCK")

#draw a picture
pdf(file = "cluster-specific-TF-heatmap.pdf",width = 8,height = 14)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[top3tfgene,], name="Regulon activity", row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled[top3tfgene,])[row_order(hm)] # to save the clustered regulons for later
dev.off()

top3gene_cell_exp <- AverageExpression(sce,
                                       assays = 'RNA',
                                       features = top3tfgene,
                                       group.by = 'seurat_clusters',
                                       slot = 'data') 
top3gene_cell_exp <- as.data.frame(top3gene_cell_exp$RNA)
top3marker_exp <- t(scale(t(top3gene_cell_exp),scale = T,center = T))
pdf(file = "TF_heatmap-expression_cluster1.pdf",width = 8,height = 8)
hm=pheatmap(top3marker_exp,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
            cluster_rows = T,cluster_cols = F,scale = "row",
            annotation_color = col,
            legendName="Relative value")
dev.off()



