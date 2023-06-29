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
