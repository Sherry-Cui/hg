######## differentiation states among cells of epithelial, neural, and mesenchymal lineages
######## CytoTRACE  ---------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggsci)
library(CytoTRACE)

load(file = './rdata/cca.sample.integrated.umap.RData')
col=pal_igv('default',alpha = 1)(51)
# mesenchymal lineage
sce <- sample.integrated[, sample.integrated$cell_type %in% c('hPSC','Mesenchymal#1','Mesenchymal#2','Mesenchymal#3')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 1000)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = pheno,emb = VEC)

VlnPlot(sce, features = c('COL1A2', 'COL3A1', 'SNAI1', 'SNAI2', 'TWIST1'),
        pt.size = 0,assay = 'RNA',cols = col,flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()
VlnPlot(sce, features = c('COL6A1','FN1','FBLN1','LAMA4','FBN2','LAMB1',
                                  'COL1A1','EFEMP1','COL1A2','BGN','COL2A1','COL3A1'),
        pt.size = 0,assay = 'RNA',cols = col,flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()
VlnPlot(sce, features = c('RSPO3','WNT5A','CXCL12'),pt.size = 0,assay = 'RNA',
        cols = col,flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()
saveRDS(sce,file = 'mesen.rds')

# neural lineage
sce <- sample.integrated[, sample.integrated$cell_type %in% c('hPSC','Early NPC','Premigratory ENC','NPC','Migratory ENCC','Neuron')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 1000)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = pheno,emb = VEC)

VlnPlot(sce, features = c('SOX2','SOX11', 'PAX6', 'SOX10', 'CDH6', 'NGFR', 
                             'MAP2'),pt.size = 0,assay = 'RNA',cols = col,
        flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()
VlnPlot(sce, features = c('COL1A2','FN1','LAMB1','FBLN2','TNC','COL1A1','COL5A2',
                             'VCAN','FBLN1','LAMA5','FBN2','COL2A1','LAMA4',
                             'COL11A2','POSTN','LGI4'),pt.size = 0,assay = 'RNA',
        cols = col,flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()
VlnPlot(sce, features = c('WNT5A','TGFB2','DLL1','DLL2',
                             'DLL3','DKK3'),pt.size = 0,assay = 'RNA',cols = col,
        flip = T,stack = T,group.by = 'cell_type',split.by = 'cell_type')+NoLegend()

saveRDS(sce,file = 'neuron.rds')

# epithelial lineage
sce <- sample.integrated[, sample.integrated$cell.type %in% c('hPSC','DE','Partial Epi','Epithelium','Enteroendocrine')]
DefaultAssay(sce) <- 'integrated'
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = 'umap')

data <- as.matrix(GetAssayData(sce, assay = "RNA", slot = "counts"))
results <- CytoTRACE(data, ncores = 12, subsamplesize = 500)
pheno <- sce$lab
pheno <- as.character(pheno)
names(pheno) <- rownames(sce@meta.data)
VEC=sce@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = pheno,emb = VEC)

VlnPlot(sce, features = c('COL1A2','FBLN1','COL4A6','COL4A2','LAMC1','MATN2',
                          'FBLN2','COL12A1','FN1','COL18A1','LAMB1','HSPG2','COL2A1','VTN','VCAN'),pt.size = 0,assay = 'RNA',cols = col,ncol = 5)
VlnPlot(sce, features = c('POU5F1','SOX17', 'GATA4', 'SOX2', 'NKX2-2'),pt.size = 0,
        assay = 'RNA',cols = col,flip = T,stack = T,group.by = 'cell.type',split.by = 'cell.type')+NoLegend()
VlnPlot(sce, features = c('COL1A2','FBLN1','COL4A6','COL4A2','LAMC1','MATN2',
                          'FBLN2','FN1','COL18A1','LAMB1','HSPG2','COL2A1','VTN','VCAN'),
        pt.size = 0,assay = 'RNA',cols = col,flip = T,stack = T,group.by = 'cell.type',split.by = 'cell.type')+NoLegend()

saveRDS(sce,file = 'epi.rds')


############# Pseudotime and pseudospace of Epithelium cells 
library(monocle)
library(URD)
library(clusterProfiler)
library(org.Hs.eg.db)
###### monocle2 ---------------------------------------------------------------
load(file = 'epi.RData')
DefaultAssay(counts) <- 'RNA'
data <- as(as.matrix(counts@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = counts@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

######### GO analysis ---------------------------------------------------------
# antrum ---------
load(file = 'scanpy_reduction_antral.rdata')
mycds_antrum <- mycds[,pData(mycds)$lab %in% "Antral Epi"]
pData(mycds_antrum)$Pseudotime <- antral$dpt_pseudotime

pseudotime_de <- differentialGeneTest(mycds_antrum, fullModelFormulaStr = "~sm.ns(Pseudotime)") 
pseudotime_de <- subset(pseudotime_de, qval < 1e-7) 
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval),]
p = plot_pseudotime_heatmap(mycds_antrum[pseudotime_de$gene_short_name,],return_heatmap=T,
                            hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),show_rownames = F,cluster_rows = T,num_clusters = 4)
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- pseudotime_de[hp.genes, c("gene_short_name", "pval", "qval")]
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
colnames(clustering) <- "Gene_Clusters"
table(clustering)
clustering$gene_short_name <- rownames(clustering)
mg <- clustering[hp.genes,c('Gene_Clusters','gene_short_name')]
mg <- cbind(mg,Time_diff_sig)
mg <- mg[,-2]
mg$Gene_Clusters <- plyr::mapvalues(x=mg$Gene_Clusters, from=c('1','2','4','3'),
                                    to=c('1','2','3','4'))
mg <- mg[order(mg$Gene_Clusters),]
p = plot_pseudotime_heatmap(mycds_antrum[mg$gene_short_name,],return_heatmap=T,
                            hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),show_rownames = F,cluster_rows = F)

ids=bitr(mg$gene_short_name,'SYMBOL','ENTREZID','org.Hs.eg.db') 
markers=merge(mg,ids,by.x='gene_short_name',by.y='SYMBOL') 
markers <- markers[!duplicated(markers$gene_short_name),]
rownames(markers) <- markers$gene
markers$module <- markers$Gene_Clusters
markers$module <- plyr::mapvalues(x=markers$module, from=c("1","2","3",'4'), 
                                  to=c("Module1","Module2","Module3",'Module4'))
gcSample=split(markers$ENTREZID,markers$module)
bp <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
dotplot(bp) + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
))
# fundus ---------
load(file = 'scanpy_reduction_fundus.rdata')
mycds_fundus <- mycds[,pData(mycds)$lab %in% "Fundic Epi"]
pData(mycds_fundus)$Pseudotime <- fundus$dpt_pseudotime
fundus_pseudotime_de <- differentialGeneTest(mycds_fundus, fullModelFormulaStr = "~sm.ns(Pseudotime)") 
fundus_pseudotime_de <- fundus_pseudotime_de[order(fundus_pseudotime_de$qval),]
fundus_pseudotime_de <- subset(fundus_pseudotime_de, qval < 1e-4) 
p = plot_pseudotime_heatmap(mycds_fundus[fundus_pseudotime_de$gene_short_name,],return_heatmap=T,
                            hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),show_rownames = F,cluster_rows = T,num_clusters = 4)
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- fundus_pseudotime_de[hp.genes, c("gene_short_name", "pval", "qval")]
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
colnames(clustering) <- "Gene_Clusters"
table(clustering)
clustering$gene_short_name <- rownames(clustering)
mg <- clustering[hp.genes,c('Gene_Clusters','gene_short_name')]
mg <- cbind(mg,Time_diff_sig)
mg <- mg[,-2]
mg$Gene_Clusters <- plyr::mapvalues(x=mg$Gene_Clusters, from=c('3','1','4','2'),
                                    to=c('1','2','3','4'))
mg <- mg[order(mg$Gene_Clusters),]
p = plot_pseudotime_heatmap(mycds_fundus[mg$gene_short_name,],return_heatmap=T,
                            hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),show_rownames = F,cluster_rows = F)

ids=bitr(mg$gene_short_name,'SYMBOL','ENTREZID','org.Hs.eg.db') 
markers=merge(mg,ids,by.x='gene_short_name',by.y='SYMBOL') #
markers <- markers[!duplicated(markers$gene_short_name),]
rownames(markers) <- markers$gene
markers$module <- markers$Gene_Clusters
markers$module <- plyr::mapvalues(x=markers$module, from=c("1","2","3",'4'), 
                                  to=c("Module1","Module2","Module3",'Module4'))
gcSample=split(markers$ENTREZID,markers$module)
bp <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
dotplot(bp) + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
))

######## URD on D16 epithelial cells -----------------------------------------
load(file = 'epiremovegland.rdata')
d16 <- epiremovegland[, epiremovegland$day %in% "D16" ] #12447
d16$cell.type <- ordered(d16$cell.type,levels = c("Fundus#1","Fundus#2","Antrum#1","Antrum#2","Antrum#3"))

de.gland <- d16
count <- de.gland@assays$RNA@counts
meta <- de.gland@meta.data
object <- createURD(count.data=count, meta=meta, min.cells = 20, min.counts=20)
object <- createURD(count.data=count, meta=meta, min.cells = 0, min.counts=0,min.genes = 0,gene.max.cut = 15000,max.genes.in.ram = 15000)
rm(list=c("count", "meta"))
shhhh <- gc()
stages <- unique(object@meta$cell.type)
cells.each.stage <- lapply(stages, function(stage) rownames(object@meta)[which(object@meta$cell.type == stage)])
var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(object, cells.fit=cells.each.stage[[n]], set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=stages[[n]], do.plot=T))
names(var.genes.by.stage) <- stages
var.genes <- sort(unique(unlist(var.genes.by.stage))) 
object@var.genes <- var.genes

object <- calcPCA(object)
set.seed(18)
object <- calcTsne(object, perplexity = 30, theta=0.5)
set.seed(17)
object <- graphClustering(object, dim.use="pca", num.nn=c(15,20,30), do.jaccard=T, method="Louvain")
object <- calcKNN(object, nn=100)

object <- calcDM(object, knn=200, sigma.use=8)
root.cells <- rownames(object@meta)[object@meta$cell.type=="Fundus#1"]
flood.result <- floodPseudotime(object, root.cells=root.cells, n=10, minimum.cells.flooded=2, verbose=T)
object <- floodPseudotimeProcess(object, flood.result, floods.name="pseudotime", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)
gg.data <- cbind(object@pseudotime, object@meta[rownames(object@pseudotime),])
ggplot(gg.data, aes(x=pseudotime, color=cell.type, fill=cell.type)) + geom_density(alpha=0.4) + theme_bw()

cds_sub <- mycds[,pData(mycds)$cell.type %in% c("Antral Epi","Fundic Epi")]
cds_sub <- cds_sub[,pData(cds_sub)$day %in% 'D16']
pse <- as.data.frame(object@pseudotime)
pse$cell <- rownames(pse)
pse[is.na(pse)] <- 0
pdata <- pData(cds_sub)
pdata$cell <- rownames(pdata)
mg <- full_join(pdata,pse,by = 'cell')
pData(cds_sub)$Pseudotime <- mg$pseudotime
gene <- c('HOXA-AS3','HOXB6','HOXD-AS2','HOXA7','HOXB-AS3','HOXD3','HOXD1','HOXD4','HOXB1','HOXB9','HOXB-AS1','HOXC6',
          'HOXD9','HOXA1','HOXB8','HOXA5','HOXC4','HOXA4','HOXB5','HOXB4',
          'HOXA2','HOXB2','HOXA3','HOXB3','HOXB7','HOXA-AS2','HOXC5')
pseudotime_de <- differentialGeneTest(cds_sub[gene,], fullModelFormulaStr = "~sm.ns(Pseudotime)") 
plot_pseudotime_heatmap(cds_sub[pseudotime_de$gene_short_name,],return_heatmap=T, hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),show_rownames = T,cluster_rows = F)












