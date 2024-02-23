######## compare with the data in vivo 
library(Matrix)
library(Seurat)
library(ggsci)
library(pheatmap)

invivo <- readRDS(file = 'Dat_HIO_Cl12_and_d132_stem_cell_fetal_pt_gene_average_expr_as_ref_for_qp.rds')
index <- read.csv(file = 'Table_fetal_atlas_cell_index.csv',header = F)
symbol <- read.csv(file = 'Table_fetal_atlas_gene_symbol.csv',header = F)
info <- read.csv(file = 'Table_fetal_atlas_meta_info.csv')
count <- Matrix::readMM('Table_fetal_atlas_count.mtx')
load(file = 'sample.integrated.RData')

rownames(count) <- symbol$V1 
colnames(count) <- index$V1
all <- CreateSeuratObject(count,meta.data = info) #28691 155232
all <- all %>%
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(dims = 1:20)%>%
  RunUMAP(dims = 1:20)
cell.embeddings <- info[,c('UMAP_X','UMAP_Y')]
cell.embeddings <- as.matrix(cell.embeddings)
colnames(cell.embeddings)[1] <- 'UMAP_1'
colnames(cell.embeddings)[2] <- 'UMAP_2'
all@reductions$umap@cell.embeddings <- cell.embeddings
col <- levels(as.factor(all$cluster_col))
DimPlot(all, reduction = "umap", group.by = "Corrected_tissue",raster=FALSE,split.by = 'Corrected_tissue',ncol = 4,label = F,cols = col) 

stomach <- all[, all$Corrected_organ_group %in% 'Stomach']
p1 <- DimPlot(stomach, reduction = "umap", group.by = "Age_week",raster=FALSE,label = F,cols = col[15:27])
p2 <- DimPlot(stomach, reduction = "umap", group.by = "Major_cell_type",raster=FALSE,label = TRUE,cols = col) 
p3 <- DimPlot(stomach, reduction = "umap", group.by = "Cell_type",raster=FALSE,label = F,cols = col)
p1|p2|p3

# addmodulescore
Idents(stomach) <- 'Major_cell_type'
markers <- FindAllMarkers(stomach, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Mesenchymal <- list(subset(markers,cluster=='Mesenchymal')$gene)
DefaultAssay(sample.integrated) <- 'RNA'
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = Mesenchymal,
  name = 'lite_Mesen'
)
gs=list(
  Epithelial = subset(markers,cluster=='Epithelial')$gene, 
  Endothelial = subset(markers,cluster=='Endothelial')$gene, 
  Neuronal = subset(markers,cluster=='Neuronal')$gene
)
gs = lapply(gs, toupper)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = gs
)
p1 <- FeaturePlot(sample.integrated,features = 'lite_Mesen1',cols = c('lightgrey','red'),label = F) +ggtitle('lite_Mesenchymal')
p2 <- FeaturePlot(sample.integrated,features = 'Cluster1',cols = c('lightgrey','red'),label = F) +ggtitle('lite_Epithelial')
p3 <- FeaturePlot(sample.integrated,features = 'Cluster2',cols = c('lightgrey','red'),label = F) +ggtitle('lite_Endothelial')
p4 <- FeaturePlot(sample.integrated,features = 'Cluster3',cols = c('lightgrey','red'),label = F) +ggtitle('lite_Neuronal')
(p1|p2)/(p3|p4)

# heatmap 
set.seed(20230518)
sample <- subset(sample.integrated,downsample=2000) #23068
list <- list(stomach,sample)
for (i in 1:length(list)){
  list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
  list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", nfeatures = 2000, verbose = T)
}
anchors <- FindIntegrationAnchors(object.list = list,k.filter=200, anchor.features = 2000)
integrated <- IntegrateData(anchorset = anchors)
integrated <- ScaleData(integrated, verbose = T)
integrated <- RunPCA(integrated,verbose = T)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated$lab1[1:13372] <- integrated$Major_cell_type[1:13372]
integrated$lab1[1:13372] <- paste('lite',integrated$lab1[1:13372],sep = '_')
integrated$lite <- 'in vitro'
integrated$lite[1:13372] <- 'in vivo'

cols <- c("#ff4a46","#008941","#ffdbe5","#0000a6","#eec3ff","#456d75")
Idents(vitro) <- 'day'
d16 <- vitro[, Idents(vitro) %in% "D16"] 
p1 <- DimPlot(d16, group.by = "lab1", pt.size=0.1,label = T,reduction = 'umap')+
  xlim(-13, 10) + ggtitle("In vitro(D16)") + scale_color_manual(values = col)
p2 <- DimPlot(vivo, group.by = "lab1", pt.size=0.1,label = T) + ggtitle("In vivo") + scale_color_manual(values = cols)
p1|p2

DefaultAssay(d16) <- 'RNA'
DefaultAssay(vivo) <- 'RNA'
p1 <- FeaturePlot(d16,features = c('CDH1','COL1A2','MAP2','FLT1'),cols = c('lightgrey','red'),ncol = 1,label = F)
p1 <- FeaturePlot(d16,features = 'CDH1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p2 <- FeaturePlot(d16,features = 'COL1A2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p3 <- FeaturePlot(d16,features = 'MAP2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p4 <- FeaturePlot(d16,features = 'FLT1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p <- p1/p2/p3/p4
p5 <- FeaturePlot(vivo,features = c('CDH1','COL1A2','MAP2','FLT1'),cols = c('lightgrey','red'),ncol = 1,label = F)
p5 <- FeaturePlot(vivo,features = 'CDH1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p6 <- FeaturePlot(vivo,features = 'COL1A2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p7 <- FeaturePlot(vivo,features = 'MAP2',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p8 <- FeaturePlot(vivo,features = 'FLT1',cols = c('lightgrey','red'),ncol = 1,label = F)+ xlim(-13, 10)
p|(p5/p6/p7/p8)









