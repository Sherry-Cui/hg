## Integrate all 3D samples 
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)

load(file = "D4.O.S1.1.filter.rdata") #10925
load(file = "D4.O.S1.2.filter.rdata") #10723
load(file = "D4.O.S2.1.filter.rdata") #18813
load(file = "D4.O.S2.2.filter.rdata") #16460

load(file = "D7.O.S1.1.filter.rdata") #10488
load(file = "D7.O.S1.2.filter.rdata") #9154
load(file = "D7.O.S2.1.filter.rdata") #11206
load(file = "D7.O.S2.2.filter.rdata") #10398


load(file = "D10.O.S1.1.filter.rdata") #5480
load(file = "D10.O.S2.1.filter.rdata") #4545
load(file = "D10.O.S2.2.filter.rdata") #5154

load(file = "D13.O.S1.1.filter.rdata") #12053
load(file = "D13.O.S1.2.filter.rdata") #11975
load(file = "D13.O.S2.1.filter.rdata") #13230
load(file = "D13.O.S2.2.filter.rdata") #11605

load(file = 'PDMS_1.filter.rdata') #10926
load(file = 'PDMS_2.filter.rdata') #12935
load(file = 'dish.filter.rdata') #8484


sample.list <- list(D4.O.S1.1.f,D4.O.S1.2.f,D4.O.S2.1.f,D4.O.S2.2.f,
                    D7.O.S1.1.f,D7.O.S1.2.f,D7.O.S2.1.f,D7.O.S2.2.f,D10.O.S1.1.f,D10.O.S2.1.f,
                    D10.O.S2.2.f,D13.O.S1.1.f,D13.O.S1.2.f,D13.O.S2.1.f,D13.O.S2.2.f,PDMS_1.f,PDMS_2.f,dish.f)
for (i in 1:length(sample.list)){
  sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
  sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- ScaleData(x,verbose = FALSE)
  x <- RunPCA(x,verbose = FALSE)
})

sample.anchors <- FindIntegrationAnchors(object.list = sample.list,k.filter=200, anchor.features = 2000, verbose = F)
sample.integrated <- IntegrateData(anchorset = sample.anchors, verbose = F) #266396

DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- ScaleData(sample.integrated, verbose = F)
sample.integrated <- RunPCA(sample.integrated,verbose = F)
ElbowPlot(sample.integrated, ndims = 50)

sample.integrated <- RunUMAP(sample.integrated,reduction = "pca", dims = 1:20, verbose = F)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:20, verbose = F)
sample.integrated <- FindClusters(sample.integrated, resolution = seq(0.1,1,0.1), verbose = F)

save(sample.integrated,file = "sample.integrated.RData")

######## Plot  -------------------------------------------------------------
col=pal_igv('default',alpha = 1)(51)
p1 <- DimPlot(sample.integrated, reduction = "umap", group.by = "day",label = F,raster=FALSE)+ scale_color_manual(values = col[35:51])
p2 <- DimPlot(sample.integrated, reduction = "umap", group.by = "cell.type",label = F,raster=FALSE)+ scale_color_manual(values = col)
p1|p2

Cellratio <- prop.table(table(sample.integrated$cell.type, sample.integrated$day), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$cell_type <- Cellratio$Var1
Cellratio$day <- Cellratio$Var2
ggplot(Cellratio) + 
  geom_bar(aes(x =day, y= Freq, fill = cell_type),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Day',y = 'Cell type')+
  scale_fill_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

gene <- c('POU5F1','NANOG', #hPSC
          'SOX17','EOMES','GSC','MIXL1',#DE
          'CLDN4','CLDN18','GATA4','CDH1',#Epithelium
          'COL1A2','COL3A1',#Mesenchymal
          'OLIG2','PCGF6',#Early NPC
          'SOX2','PAX6','NES',#NPC
          'SOX10','PAX3','EDNRB',#Premigratory ENNC
          'SOX11','CDH6','B3GAT1',#migratory ENNC
          'TUBB3','MAP2','NGFR',#Neuron
          'CDH5','FLT1', #Endothelial
          'SST','GHRL')#Enteroendocrine
DotPlot(sample.integrated, features =gene,cols = c("lightgrey",'#FF0000'))+ RotatedAxis() 



######## subset epithelial cells 
dt = sample.integrated[, sample.integrated$cell.type %in% "Epithelium" ]
counts <- dt@assays$RNA@counts
counts <- CreateSeuratObject(counts = counts)
counts$orig.ident <- dt$orig.ident
counts$day <- dt$day
list <- SplitObject(counts, split.by = "orig.ident")
for (i in 1:length(list)){
  list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
  list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = 2000)
counts <- IntegrateData(anchorset = anchors) #26186
DefaultAssay(counts) <- "integrated"
counts <- counts %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20)
counts <- FindNeighbors(counts, reduction = "pca", dims = 1:20, verbose = F)
counts <- FindClusters(counts,resolution = 0.4, verbose = F)
counts <- FindClusters(counts,resolution = 0.1, verbose = F)
counts <- FindClusters(counts,resolution = 0.3, verbose = F)
counts <- FindClusters(counts,resolution = 0.2, verbose = F)
gene <- c('PDX1','GAST','MUC5AC','CLDN18','MUC6','KLF5',# Antrum
          'ATP4B','ATP4A','PGA3','PGC','IRX2','IRX3','IRX5','BHLHA15','MUC2','GHRL',# Fundus
          'TFF2','TFF1','TFF3')# gland                                
DefaultAssay(counts) <- 'RNA'
DotPlot(counts, features =gene,cols = c("lightgrey",'#FF0000'))+ RotatedAxis() 

Cellratio <- prop.table(table(counts$cell.type, counts$day), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$cell_type <- Cellratio$Var1
Cellratio$day <- Cellratio$Var2
ggplot(Cellratio) + 
  geom_bar(aes(x =day, y= Freq, fill = cell_type),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Day',y = 'Cell type')+
  scale_fill_manual(values = col[4:51])+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

# DEG
epiremovedgland <- counts[, counts$cell.type %in% c("Fundus#1","Fundus#2","Antrum#1","Antrum#2","Antrum#3")] 
markers <- FindAllMarkers(epiremovedgland, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(epiremovedgland) <- 'RNA'
epiremovedgland <- ScaleData(epiremovedgland)
mt <- as.data.frame(t(as.matrix(GetAssayData(epiremovedgland, assay = "RNA", slot = "scale.data"))))
group.by <- 'cell.type'
mt <- aggregate(mt, by=list(epiremovedgland@meta.data[[group.by]]), FUN="mean")
rownames(mt) <- mt$Group.1
mt <- t(mt[,-1])
markers %>%
  group_by(cluster) %>%
  top_n(n = 10 ,wt = avg_log2FC) -> top10
cts <- as.matrix(mt[top10$gene,])
ComplexHeatmap::pheatmap(cts,show_colnames =T,show_rownames = T,
                         color =colorRampPalette(rev(brewer.pal(n = 35, name ="RdYlBu")))(100),
                         cluster_rows = F,
                         cluster_cols = F,
                         name= 'Scaled Expression')

save(counts,file = "epi.RData")
save(counts,file = "epiremovedgland.RData")













