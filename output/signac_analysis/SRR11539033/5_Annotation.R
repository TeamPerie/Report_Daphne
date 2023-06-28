library(Signac)
library(Seurat)
library(readr)

SO <- read_rds("SO_LabelTransfer.rds")
SO <- RenameIdents(
  object = SO,
  '0' = 'Memory-CD4-TCell',
  '1' = 'GammaDelta-TCell',
  '2' = 'Activated-BCell',
  '3' = 'Naive-CD4-TCell',
  '4' = 'NK-Cell',
  '5' = 'NK-Cell',
  '6' = 'Cytotoxic-CD8-TCell',
  '7' = 'CD14-Monocyte',
  '8' = 'NK-Cell',
  '9' = 'CD14-Monocyte',
  '10' = 'Memory-BCell',
  '11' = 'DC',
  '12' = 'Treg',
  '13' = 'CD16-Monocyte',
  '14' = 'DC'
)
SO <- AddMetaData(SO, metadata = Idents(SO), col.name = "Cell.Types")
write_rds(SO, "SO_Annotated.rds")

png("UMAP.png")
p <- DimPlot(object = SO, label = TRUE) + NoLegend()
print(p)
dev.off()

###lymphoid/myeloid
SO <- read_rds("SO_Annotated.rds")
SO <- RenameIdents(
  object = SO,
  'Memory-CD4-TCell' = 'Lymphoid',
  'GammaDelta-TCell' = 'Lymphoid',
  'CD14-Monocyte' = 'Myeloid',
  'Activated-BCell' = 'Lymphoid',
  'Naive-CD4-TCell' = 'Lymphoid',
  'NK-Cell' = 'NK-Cell',
  'Cytotoxic-CD8-TCell' = 'Lymphoid',
  'CD16-Monocyte' = 'Myeloid',
  'Memory-BCell' = 'Lymphoid',
  'DC' = 'DC',
  'Treg' = 'Lymphoid'
)

png("UMAP_Lineages.png")
DimPlot(object = SO, label = TRUE, cols = c("#F8766D", "#00BFC4", "grey", "grey")) + NoLegend()
dev.off()
