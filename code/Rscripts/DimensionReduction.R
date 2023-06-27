library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
library(SeuratData)
library(randomcoloR)
library(dplyr)
library(BuenColors)

args = commandArgs(trailingOnly = TRUE)

SO <- read_rds("SO_QCFiltered.rds")

#Normalization and dimension reduction
#Latent Semantic Indexing (LSI)
SO <- RunTFIDF(SO) #log(TFxIDF) and scale factor 10000
SO <- FindTopFeatures(SO) #based on total count, no min.cutoff
SO <- RunSVD(SO) #Use all features, scale cell embedding to mean 0 and SD 1

#Check if first dimension is related to read depth
png("depth_dimension_correlation.png")
p <- DepthCor(SO)
print(p)
dev.off()

#Clustering and UMAP, without first dimension
SO <- RunUMAP(SO, reduction = "lsi", dims = 2:30)
SO <- FindNeighbors(SO, reduction = "lsi", dims = 2:30, k.param = 15)
SO <- FindClusters(SO, resolution = 0.85, algorithm = 1) #original Louvain algorithm

png("Clustered_UMAP.png")
p <- DimPlot(object = SO, label = TRUE) + NoLegend()
print(p)
dev.off()

#Compute gene accessibility
gene.activities <- GeneActivity(SO)

# add to the Seurat object as a new assay
SO[['RNA']] <- CreateAssayObject(counts = gene.activities)
SO <- NormalizeData(SO, assay = 'RNA', normalization.method = 'LogNormalize') #scaling factor 10000
SO <- FindVariableFeatures(SO, assay="RNA") #vst method
SO <- ScaleData(SO, assay="RNA")

if (args == "BM"){
  # Import RNA data, preprocess and only keep HSC and progenitor cells (CD34+)
  # Find anchors in shared space reduced by Canonical Correlation Analsysis
  InstallData("bmcite")
  bm.rna <- LoadData(ds = "bmcite")
  Idents(bm.rna) <- bm.rna@meta.data$celltype.l2
  bm.rna <- subset(bm.rna, idents = c("Prog_RBC", "GMP", "Prog_Mk", "Prog_DC", "Prog_B 1", "pDC", "HSC", "Prog_B 2", "LMPP"))
  bm.rna <- NormalizeData(bm.rna, normalization.method = "LogNormalize", assay="RNA")
  bm.rna <- FindVariableFeatures(bm.rna, assay="RNA")
  bm.rna <- ScaleData(bm.rna, assay="RNA")
  
  # Do label transfer
  transfer.anchors <- FindTransferAnchors(reference = bm.rna, query = SO, features = VariableFeatures(object = bm.rna), 
                                          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = bm.rna$celltype.l2,
                                       weight.reduction = SO[['lsi']],
                                       dims = 2:30)
}
if (args == "PBMC"){
  #Derived from https://github.com/caleblareau/mtscATACpaper_reproducibility/blob/master/figure_paired_cd34_pbmc/code/04_PBMC_annotation.R
  pbmc.rna <- readRDS("../../../data/reference/13March2020_recluster_reannotated_10xv3.rds")
  pbmc.rna <- pbmc.rna[,pbmc.rna@meta.data$celltype != "IFNactive_monocyte"] # Remove ambiguous myeloid/Tcell cluster
  
  # Do label transfer
  transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = SO, features = VariableFeatures(object = pbmc.rna), 
                                          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype,
                                       weight.reduction = SO[['lsi']],
                                       dims = 2:30)
}
SO <- AddMetaData(SO, metadata = celltype.predictions)
write_rds(SO, "SO_LabelTransfer.rds")

#Plot predicted cell types
palette <- distinctColorPalette(length(unique(unique(SO@meta.data$predicted.id))))
png("TransferCelltype.png", width = 900, height = 600)
p <- DimPlot(
  object = SO,
  group.by = "predicted.id",
  pt.size = 1.5,
  cols = palette
)
print(p)
dev.off()

#Assignment probabilities
#The proportion of cells in each cluster assigned to each cell type
long_df <- SO@meta.data %>% group_by(seurat_clusters, predicted.id) %>%
  summarize(count = n()) %>% ungroup() %>%
  group_by(seurat_clusters) %>% mutate(prop = count / sum(count))
reshape2::dcast(long_df, seurat_clusters ~ predicted.id, value.var = "prop", fill = 0)
long_df$predicted.id <- factor(as.character(long_df$predicted.id), rev(unique(as.character(long_df$predicted.id))))

plot <- ggplot(long_df, aes(y = predicted.id, x = seurat_clusters, fill = prop)) + 
  geom_tile() + pretty_plot(fontsize = 6) + L_border() +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos")) +
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "none") + labs(x = "ATAC Cluster", y = "scRNA transfer label")
cowplot::ggsave2(plot, file = "Heatmap_Celltypes.png", width = 2.1, height = 3)
