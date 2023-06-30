

## Annotate the cells based on the given lineages
#
# df = dataframe containing cells, mutations and annotations
# lineage = the lineages to test
# remove_other = remove non-lineage cells from the analysis?
##
annotate_lineages <- function(df, lineages, remove_other = FALSE){
  if (all(lineages == c("Lymphoid", "Myeloid"))){
      df$annotation[df$annotation %in% c("CD16-Monocyte", "CD14-Monocyte")] <- "Myeloid"
      df$annotation[df$annotation %in% c("Activated-BCell", "GammaDelta-TCell", "Naive-CD4-TCell", 
                                        "Memory-BCell", "Treg", "Cytotoxic-CD8-TCell", "Memory-CD4-TCell")] <- "Lymphoid"
  }
  if (all(lineages == c("HSC", "Other"))){
      df$annotation[df$annotation != "HSC"] <- "Other"
  }        
  if (remove_other){df <- df[df$annotation %in% lineages,]}
  return(df)	
}


library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
lineages <- strsplit(args[1], "_")[[1]]
df <- fread("table_cellular_mutation.tsv", sep = "\t", header = TRUE, data.table = FALSE)
df <- annotate_lineages(df, lineages, remove_other = FALSE)

PBMC1 <- readRDS("SO_PBMC1_variants.rds")
PBMC1_UMAP <- as.data.frame(PBMC1@reductions$umap@cell.embeddings)
PBMC1_UMAP %>% mutate(cell = rownames(PBMC1_UMAP)) -> PBMC1_UMAP

df2 <- inner_join(PBMC1_UMAP, df[,c("cell","annotation")], by = "cell") 

colors <- c("Lymphoid" = "red", "Myeloid" = "blue", "NK-Cell" = "grey", "DC" = "grey")
png("UMAP_lineages.png")
p <- ggplot(df2, aes(x = UMAP_1, y = UMAP_2, color = annotation)) +
  geom_point() +
  scale_color_manual(values = colors) +
  theme_minimal()
print(p)
dev.off()