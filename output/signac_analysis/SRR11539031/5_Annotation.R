library(Signac)
library(Seurat)
library(readr)

SO <- read_rds("SO_LabelTransfer.rds")
SO <- RenameIdents(
  object = SO,
  '0' = 'GMP',
  '1' = 'LMPP',
  '2' = 'GMP',
  '3' = 'prog-RBC',
  '4' = 'HSC',
  '5' = 'prog-RBC',
  '6' = 'prog-RBC',
  '7' = 'prog-B',
  '8' = 'prog-B',
  '9' = 'pDC',
  '10' = 'prog-B',
  '11' = 'GMP',
  '12' = 'prog-MK'
)
SO <- AddMetaData(SO, metadata = Idents(SO), col.name = "Cell.Types")
write_rds(SO, "SO_Annotated.rds")


png("UMAP.png")
p <- DimPlot(object = SO, label = TRUE) + NoLegend()
print(p)
dev.off()