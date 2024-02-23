######## CellChat of day 10 cells -------------------------------------------
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)

load(file = 'sample.integrated.RData')
sce <- sample.integrated[, sample.integrated$day %in% 'D10'] 
sce <- sce[, sce$cell.type %in% c("Fundic Epi",'Antral Epi',"Mesenchymal","Early NPC","NPC","Neuron", "ENCC")]
data.input  <- sce@assays$RNA@data
identity = data.frame(group =sce$cell.type , row.names = names(sce$cell.type)) 
unique(identity$group) 
cellchat <- createCellChat(data <- data.input) 
summary(cellchat)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "cell.type")
cellchat <- setIdent(cellchat, ident.use = "cell.type") 
levels(cellchat@idents) 
table(cellchat@idents)

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)

netVisual_aggregate(cellchat, signaling = 'MK',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'PTN',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'MIF',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'ncWNT',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'GDF',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'IGF',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'PARs',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'SPP1',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'SEMA3',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = 'FGF',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'BMP',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'WNT', vertex.receiver = vertex.receiver, pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'GRN',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'SPP1',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'CXCL',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'GAS',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'PDGF',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'VISFATIN',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = 'KIT',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'NPR2',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'NT',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'ENHO',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = 'PROS',  vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 1.7)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# figure6
netVisual_aggregate(cellchat, signaling = 'WNT',layout = 'hierarchy',pt.title=2,vertex.label.cex = 0.5)
netVisual_heatmap(cellchat, signaling = 'WNT', color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = 'WNT')
netVisual_bubble(cellchat, remove.isolate = FALSE, signaling = "WNT")

netAnalysis_signalingRole_network(cellchat, signaling = 'PDGF', width = 16, height = 8, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = 'BMP', width = 16, height = 8, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = 'WNT', width = 16, height = 8, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = 'ncWNT', width = 16, height = 8, font.size = 10)


netVisual_bubble(cellchat, remove.isolate = FALSE, signaling = c("PDGF",'BMP'), 
                 sources.use = c('Antral Epi',"Fundic Epi"),targets.use =c('Mesenchymal','Early NPC','NPC','Neuron','ENCC'))

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.heatmap = "GnBu", width = 7, height = 8)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "GnBu", width = 7, height = 8)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 1), xpd = FALSE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,title.name = "Interaction weights/strength")

saveRDS(cellchat,file = 'cellchat.rds')





