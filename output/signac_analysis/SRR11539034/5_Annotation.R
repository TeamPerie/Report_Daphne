library(Signac)
library(Seurat)
library(readr)

SO <- read_rds("SO_LabelTransfer.rds")
SO <- RenameIdents(
  object = SO,
  '0' = 'NK-Cell',
  '1' = 'GammaDelta-TCell',
  '2' = 'Naive-CD4-TCell',
  '3' = 'Activated-BCell',
  '4' = 'NK-Cell',
  '5' = 'Treg',
  '6' = 'Memory-CD4-TCell',
  '7' = 'CD14-Monocyte',
  '8' = 'Cytotoxic-CD8-TCell',
  '9' = 'Memory-BCell',
  '10' = 'CD14-Monocyte',
  '11' = 'DC',
  '12' = 'CD16-Monocyte',
  '13' = 'NK-Cell'
)
SO <- AddMetaData(SO, metadata = Idents(SO), col.name = "Cell.Types")
write_rds(SO, "SO_Annotated.rds")

png("UMAP.png")
p <- DimPlot(object = SO, label = TRUE) + NoLegend()
print(p)
dev.off()