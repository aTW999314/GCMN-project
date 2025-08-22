
####cell interaction analysis####
library(CellChat)
library(ggalluvial)
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
setwd('/home/Tanwei/GCMN')
load("GCMN_NS.RData")


sce$seurat_clusters[sce$seurat_clusters == "Schwann cells/melanocytes"] <- "melanocytes"

sce$celltype[sce$celltype == "Schwann cells/melanocytes"] <- "melanocytes"
sce.mergeTEN=sce

data.input <- GetAssayData(sce.mergeTEN, assay = "RNA", slot = "data")
meta = sce.mergeTEN@meta.data# a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$Group == "GCMN"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$celltype) # check the cell labels

cellchat.GCMN <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat.GCMN  <- addMeta(cellchat.GCMN, meta = meta)
cellchat  <- setIdent(cellchat.GCMN , ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat.GCMN@idents)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat.GCMN@DB <- CellChatDB.use

options(future.globals.maxSize = 2000 * 1024^2)  
cellchat=cellchat.GCMN
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

options(future.globals.maxSize = 5 * 1024^3)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.GCMN=cellchat
save(cellchat.GCMN,file="cellchat.GCMN.RData")

setwd("/home/Tanwei/GCMN/cell-chat-before")
print(cellchat.GCMN@net$count) 

#cellchat.NS
data.input <- GetAssayData(sce.mergeTEN, assay = "RNA", slot = "data")
meta = sce.mergeTEN@meta.data
cell.use = rownames(meta)[meta$Group == "NS"] 
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
cellchat.NS <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat.NS <- addMeta(cellchat.NS, meta = meta)
cellchat <- setIdent(cellchat.NS, ident.use = "celltype") # set "labels" as default cell identity

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.human 

cellchat.NS@DB <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat=cellchat.NS
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat.NS=cellchat


pdf(file = "GCMN_NS_cellinteraction1.pdf",width = 10,height = 8)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

all_celltypes <- unique(c(
  as.character(object.list[[1]]@idents), 
  as.character(object.list[[2]]@idents)
))
unified_levels <- levels(object.list[[1]]@idents)  

for (i in 1:length(object.list)) {
  object.list[[i]]@idents <- factor(
    object.list[[i]]@idents, 
    levels = unified_levels
  )
}
unified_levels <- levels(object.list[["NS"]]@idents) 


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file = "GCMN_NS_cellinteraction_circle_fixed.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


num.link <- rowSums(cellchat.GCMN@net$count) + colSums(cellchat.GCMN@net$count) - diag(cellchat.GCMN@net$count)
weight.MinMax <- c(min(num.link), max(num.link))  


cellchat.GCMN <- netAnalysis_computeCentrality(
  cellchat.GCMN, 
  slot.name = "netP" 
)


gg <- netAnalysis_signalingRole_scatter(
  cellchat.GCMN, 
  title = "cellchat.GCMN",  
  weight.MinMax = weight.MinMax
)

pdf(file = "outgoing_incoming_GCMN.pdf",width = 5,height = 8)
patchwork::wrap_plots(gg)
dev.off()

num.link <- rowSums(cellchat.NS@net$count) + colSums(cellchat.NS@net$count) - diag(cellchat.NS@net$count)
weight.MinMax <- c(min(num.link), max(num.link))  

cellchat.NS <- netAnalysis_computeCentrality(
  cellchat.NS, 
  slot.name = "netP"  
)

gg1 <- netAnalysis_signalingRole_scatter(
  cellchat.NS, 
  title = "cellchat.NS",  
  weight.MinMax = weight.MinMax
)

pdf(file = "outgoing_incoming_NS.pdf",width = 5,height = 8)
patchwork::wrap_plots(gg1)
dev.off()


cellchat.GCMN <- netAnalysis_computeCentrality(cellchat.GCMN, slot.name = "netP")
pdf(file = "outgoing_incoming_GCMN.pdf",width = 5,height = 8)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# the slot 'netP' means the inferred intercellular communication network of signaling path
gg1 <- netAnalysis_signalingRole_scatter(cellchat.GCMN)
gg1
dev.off()

cellchat.NS <- netAnalysis_computeCentrality(cellchat.NS, slot.name = "netP")
pdf(file = "outgoing_incoming_NS.pdf",width = 5,height = 8)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# the slot 'netP' means the inferred intercellular communication network of signaling path
gg1 <- netAnalysis_signalingRole_scatter(cellchat.NS)
gg1
dev.off()

object.list <- list(GCMN = cellchat.GCMN, NS = cellchat.NS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),merge.data=TRUE)

pdf(file = "GCMN vs NS.pdf",width = 10,height = 8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
object.list <- list(GCMN = cellchat.GCMN, NS = cellchat.NS)

pdf(file = "heatmap of GCMN vs NS.pdf",width = 10,height = 8)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(file = "heatmap of GCMN vs NS_incoming.pdf",width = 10,height = 8)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(file = "heatmap of GCMN vs NS_all.pdf",width = 10,height = 8)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

levels(cellchat@idents$NS)  

new_levels <- c(
  "melanocytes",        # 1:  "melanocytes"
  "Keratinocyte",       # 2:  "Keratinocyte"
  "Fibroblast",         # 3:  "Fibroblast"
  "Endothelial cells",  # 4:  "Endothelial cells"
  "T/DC",               # 5:  "T/DC"
  "Myeloid cells",      # 6:  "Myeloid cells"
  "Mast cells"          # 7:  "Mast cells"
)


cellchat@idents$GCMN <- factor(cellchat@idents$GCMN, levels = new_levels)
cellchat@idents$NS <- factor(cellchat@idents$NS, levels = new_levels)

levels(cellchat@idents$GCMN)
levels(cellchat@idents$NS)
####Identify upregulated and downregulated signal ligand pairs####
pdf(file = "MC of GCMN vs NS_all.pdf",width = 10,height = 8)
#1=melanocyte;2=Keratinocyte;3=Fibroblast;4=Endothelial cell;5=T/DC;6=myeloid;7=mast cell
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

pdf(file = 'fb of GCMN vs NS_all.pdf',width = 10,height = 8)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

pathways.show <- c("PTN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(file = "PTN pathway of GCMN.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("IGF") 
library(ComplexHeatmap)
pdf(file = "heatmap of IGF pathway of GCMN.pdf",width = 5,height = 8)
#> Do heatmap based on a single object
par(mfrow = c(1,2), xpd=TRUE)
ht1 <- netVisual_heatmap(object.list[[1]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(object.list)[1]))
ComplexHeatmap::draw(ht1, ht_gap = unit(0.5, "cm"))
dev.off()

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NS", "GCMN")) # set factor level
pdf(file = "IGF genes of GCMN.pdf",width = 10,height = 8)
plotGeneExpression(cellchat, signaling = "IGF", split.by = "datasets", colors.ggplot = T)
dev.off()

####downsample of melanocyte in GCMN####
setwd('/home/Tanwei/GCMN')
load("GCMN_NS.RData")
sc=sce
Idents(sc) <- "celltype"
melanocytes_gcmn <- subset(sc, celltype == "melanocytes" & Group == "GCMN")
melanocytes_ns <- subset(sc, celltype == "melanocytes" & Group == "NS")

# 2. Downsample GCMN melanocytes to the number in the NS group (1005 cells)
set.seed(123)
melanocytes_gcmn_down <- melanocytes_gcmn[, sample(colnames(melanocytes_gcmn), 1005, replace = FALSE)]

# 3. Refactor GCMN samples (replace the original melanocytes with downsampled data)
other_cells_gcmn <- subset(sc, Group == "GCMN" & celltype != "melanocytes")
sc_gcmn_down <- merge(melanocytes_gcmn_down, other_cells_gcmn)
sc_ns <- subset(sc, Group == "NS") # The NS group remains unchanged
table(sc_gcmn_down@meta.data$celltype)
table(sc_ns@meta.data$celltype)

data.input <- GetAssayData(sc_gcmn_down, assay = "RNA", slot = "data")
meta = sc_gcmn_down@meta.data# a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$Group == "GCMN"]

cellchat.GCMN.down <- createCellChat(object = sc_gcmn_down, meta = meta, group.by = "celltype")
cellchat.GCMN.down  <- addMeta(cellchat.GCMN.down, meta = meta)
cellchat.GCMN.down <- setIdent(cellchat.GCMN.down, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat.GCMN.down@idents)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat.GCMN.down@DB <- CellChatDB.use
future::plan("multisession", workers = 4)

cellchat.GCMN.down <- subsetData(cellchat.GCMN.down) 
cellchat.GCMN.down <- identifyOverExpressedGenes(cellchat.GCMN.down)
cellchat.GCMN.down <- identifyOverExpressedInteractions(cellchat.GCMN.down)
cellchat.GCMN.down <- projectData(cellchat.GCMN.down, PPI.human)

options(future.globals.maxSize = 5 * 1024^3)
cellchat.GCMN.down <- computeCommunProb(cellchat.GCMN.down, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.GCMN.down <- filterCommunication(cellchat.GCMN.down, min.cells = 10)
cellchat.GCMN.down <- computeCommunProbPathway(cellchat.GCMN.down)
cellchat.GCMN.down <- aggregateNet(cellchat.GCMN.down)

save(cellchat.GCMN.down,file="cellchat.GCMN.down.RData")

setwd("/home/Tanwei/GCMN/cell-chat-before")
print(cellchat.GCMN@net$count) 

#cellchat.NS
setwd("/home/Tanwei/GCMN/cell-chat-before")
load("cellchat.NS.RData")

setwd("/home/Tanwei/GCMN/cellinteraction-downsample")
object.list_down <- list(GCMN_down = cellchat.GCMN.down, NS = cellchat.NS)
cellchat_down <- mergeCellChat(object.list_down, add.names = names(object.list_down))

# Compare the number and intensity of interactions
pdf("GCMN_down_NS_interaction.pdf", 10, 8)
gg1 <- compareInteractions(cellchat_down, group = c(1,2))
gg2 <- compareInteractions(cellchat_down, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file = "GCMN_down_NS_cellinteraction_circle_fixed.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

num.link <- rowSums(cellchat.GCMN.down@net$count) + colSums(cellchat.GCMN.down@net$count) - diag(cellchat.GCMN.down@net$count)
weight.MinMax <- c(min(num.link), max(num.link)) 

cellchat.GCMN.down <- netAnalysis_computeCentrality(
  cellchat.GCMN.down, 
  slot.name = "netP" 
)


gg <- netAnalysis_signalingRole_scatter(
  cellchat.GCMN.down, 
  title = "cellchat.GCMN",  
  weight.MinMax = weight.MinMax
)

pdf(file = "outgoing_incoming_GCMN_down.pdf",width = 5,height = 8)
patchwork::wrap_plots(gg)
dev.off()

num.link <- rowSums(cellchat.NS@net$count) + colSums(cellchat.NS@net$count) - diag(cellchat.NS@net$count)
weight.MinMax <- c(min(num.link), max(num.link))  

cellchat.NS <- netAnalysis_computeCentrality(
  cellchat.NS, 
  slot.name = "netP"  
)


gg1 <- netAnalysis_signalingRole_scatter(
  cellchat.NS, 
  title = "cellchat.NS",  
  weight.MinMax = weight.MinMax
)


pdf(file = "outgoing_incoming_NS.pdf",width = 5,height = 8)
patchwork::wrap_plots(gg1)
dev.off()


cellchat.GCMN.down <- netAnalysis_computeCentrality(cellchat.GCMN.down, slot.name = "netP")
pdf(file = "outgoing_incoming_GCMN.pdf",width = 5,height = 8)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# the slot 'netP' means the inferred intercellular communication network of signaling path
gg1 <- netAnalysis_signalingRole_scatter(cellchat.GCMN.down)
gg1
dev.off()

cellchat.NS <- netAnalysis_computeCentrality(cellchat.NS, slot.name = "netP")
pdf(file = "outgoing_incoming_NS.pdf",width = 5,height = 8)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# the slot 'netP' means the inferred intercellular communication network of signaling path
gg1 <- netAnalysis_signalingRole_scatter(cellchat.NS)
gg1
dev.off()

object.list <- list(GCMN = cellchat.GCMN, NS = cellchat.NS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),merge.data=TRUE)

pdf(file = "GCMN vs NS.pdf",width = 10,height = 8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
object.list <- list(GCMN = cellchat.GCMN.down, NS = cellchat.NS)

pdf(file = "heatmap of GCMN vs NS.pdf",width = 10,height = 8)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(file = "heatmap of GCMN vs NS_incoming.pdf",width = 10,height = 8)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(file = "heatmap of GCMN vs NS_all.pdf",width = 10,height = 8)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

levels(cellchat@idents$NS)  

new_levels <- c(
  "melanocytes",        # 1:  "melanocytes"
  "Keratinocyte",       # 2:  "Keratinocyte"
  "Fibroblast",         # 3:  "Fibroblast"
  "Endothelial cells",  # 4:  "Endothelial cells"
  "T/DC",               # 5:  "T/DC"
  "Myeloid cells",      # 6:  "Myeloid cells"
  "Mast cells"          # 7:  "Mast cells"
)

# Reset the factor levels of the GCMN grouping
cellchat@idents$GCMN <- factor(cellchat@idents$GCMN, levels = new_levels)

# Reset the factor levels of the NS group
cellchat@idents$NS <- factor(cellchat@idents$NS, levels = new_levels)

# Verify the modification results
levels(cellchat@idents$GCMN)
levels(cellchat@idents$NS)

pdf(file = "MC of GCMN vs NS_all.pdf",width = 10,height = 8)
#1=melanocyte;2=Keratinocyte;3=Fibroblast;4=Endothelial cell;5=T/DC;6=myeloid;7=mast cell
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

pdf(file = 'fb of GCMN vs NS_all.pdf',width = 10,height = 8)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
dev.off()

pathways.show <- c("PTN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(file = "PTN pathway of GCMN.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("IGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(file = "IGF pathway of GCMN.pdf",width = 10,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

#heatmap
pathways.show <- c("IGF") 
library(ComplexHeatmap)
pdf(file = "heatmap of IGF pathway of GCMN.pdf",width = 5,height = 8)
#> Do heatmap based on a single object
par(mfrow = c(1,2), xpd=TRUE)
ht1 <- netVisual_heatmap(object.list[[1]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(object.list)[1]))
ComplexHeatmap::draw(ht1, ht_gap = unit(0.5, "cm"))
dev.off()

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NS", "GCMN")) # set factor level
pdf(file = "PTN genes of GCMN.pdf",width = 10,height = 8)
plotGeneExpression(cellchat, signaling = "PTN", split.by = "datasets", colors.ggplot = T)
dev.off()
