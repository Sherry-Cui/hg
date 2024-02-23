######## figure5 TF regulons identification of day 16 epithelial subtypes and heatmap
library(SCENIC)
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(ggsci)
library(circlize)
library(stringr)
library(writexl)

load(file = 'epi.RData')
d16 <- counts[, counts$day %in% 'D16'] 
loom <- open_loom('out_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons <- regulonsToGeneLists(regulons_incidMat) 
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)     
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])   
embeddings <- get_embeddings(loom)  
close_loom(loom)
sub_regulonAUC <- regulonAUC[,match(colnames(d16),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(d16))
cellTypes <- data.frame(row.names = colnames(d16), 
                        celltype = d16$cell.type)
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes),cellTypes[,selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
melt <- melt(regulonActivity_byGroup_Scaled)
melt %>%
  group_by(Var2) %>%
  top_n(n = 10, wt = value) -> top10
plot.data <- regulonActivity_byGroup_Scaled[top10$Var1,]
p1 <- Heatmap(
  plot.data,
  name                         = "z-score",
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = T,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

# sup_figure9 regulation networks:cytoscape 
