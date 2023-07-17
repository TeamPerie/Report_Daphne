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

SO_UMAP <- as.data.frame(SO@reductions$umap@cell.embeddings)
SO_UMAP %>% mutate(cell = rownames(SO_UMAP)) -> SO_UMAP
SO_META <- as.data.frame(SO@meta.data)
SO_META %>% mutate(cell = rownames(SO_META)) -> SO_META
df <- inner_join(SO_UMAP, SO_META[,c("cell","Cell.Types")], by = "cell") 
svg("Annotated_UMAP.svg")
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cell.Types)) +
      geom_point() +
      theme_classic()
print(p)
dev.off()

