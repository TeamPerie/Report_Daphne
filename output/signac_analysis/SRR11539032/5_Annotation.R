library(Signac)
library(Seurat)
library(readr)

SO <- read_rds("SO_LabelTransfer.rds")
SO <- RenameIdents(
  object = SO,
  '0' = 'LMPP',
  '1' = 'GMP',
  '2' = 'prog-RBC',
  '3' = 'HSC',
  '4' = 'prog-RBC',
  '5' = 'GMP',
  '6' = 'prog-MK',
  '7' = 'prog-B',
  '8' = 'GMP',
  '9' = 'prog-B',
  '10' = 'prog-MK',
  '11' = 'pDC',
  '12' = '?'
)
SO <- AddMetaData(SO, metadata = Idents(SO), col.name = "Cell.Types")
write_rds(SO, "SO_Annotated.rds")

png("UMAP.png")
p <- DimPlot(object = SO, label = TRUE) + NoLegend()
print(p)
dev.off()