annotate_clones <- function(df, lineages){
    df$bias <- rep(NA, nrow(df))
    for (x in names(lineages)){
        df$bias[df$mutation %in% lineages[[x]]] <- x
    }
   df$bias[is.na(df$bias)] <- "undefined"
   if (nrow(df[df$bias == "Lineage2",]) == 0){df$bias[df$bias=="undefined"] <- "Lineage2"}
   return(df) 
}

annotate_SO <- function(df, SO){
    cells <- colnames(SO)
    new_column <- replace(cells, cells %in% df$cell[df$bias == "Lineage1"], "Lineage1")
    new_column <- replace(new_column, cells %in% df$cell[df$bias == "Lineage2"], "Lineage2")
    new_column <- replace(new_column, !(new_column %in% c("Lineage1", "Lineage2")), "undefined")
    SO <- AddMetaData(SO, new_column, col.name = "Bias")
    return(SO)   
}

create_umap <- function(SO, output_dir, Sample){
    df <- SO@reductions$umap@cell.embeddings
    df <- cbind.data.frame(df, "Bias" = SO@meta.data$Bias)
    df <- df[df$Bias != "undefined",]
    df %>% arrange(Bias) -> df
    png(paste0(output_dir, "/UMAP_",Sample,"_Lineages.png"))
    p <- ggplot(df, aes(x = -(as.numeric(UMAP_1)), y = as.numeric(UMAP_2), color = Bias)) +
    geom_point() +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Bias") + theme_classic()
    print(p)
    dev.off()
}

density_plot <-function(SO, lineages, output_dir, Sample){
    Idents(SO) <- "Bias"
    so1 <- subset(SO,idents = c("Lineage1"))
    so2 <- subset(SO,idents = c("Lineage2"))
    p1 <- createDensityPlot(so1, bins = 50, lineages[1]) + NoLegend()
    p2 <- createDensityPlot(so2, bins = 50, lineages[2])
    png(paste0(output_dir, "/", Sample, "_Density_UMAP.png"), width=480, height=960)
    print(ggarrange(p1, p2,ncol = 1, nrow = 2))
    dev.off()
}

#From https://github.com/TeamPerie/Cosgrove-et-al-2022/blob/main/Figure1/INTEGRATION/helper_methods_refactored.R
createDensityPlot <- function(sobj,bins = 10, lineage){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "umap"))
  p <- ggplot(tmp.all, aes(x = -(UMAP_1), y = UMAP_2)) + geom_point(colour="#00000000") + 
    stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=bins) + 
    scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) +
    theme_classic() + 
    theme(legend.position="bottom") +
    ggtitle(lineage) + xlim(c(-10,10)) + ylim(c(-8,8))
  
  return(p)
}

deg_between_lineages <- function(SO, output_dir, Sample){
    DefaultAssay(SO) <- "RNA"
    Idents(SO) <- SO@meta.data$Bias
    markers <- FindMarkers(SO, ident.1="Lineage1", ident.2="Lineage2", logfc.threshold = 0)
    markers$gene <- rownames(markers)
    write.table(markers, file = paste0(output_dir,"/Markers_",Sample,".tsv"), sep = "\t")
    return(markers)
}



library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
#library(clusterProfiler)
#library(enrichplot)
library(dplyr)
library(tidyr)
library(ggpubr)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
Lineages <- strsplit(args[1], "_")[[1]]
Sample <- args[2]
print(paste0("Sample: ", Sample))
print(paste0("Lineages: ", Lineages))

output_dir = getwd()
SO <- read_rds(paste0("../SO_",Sample,"_variants.rds"))
df <- fread("../table_cellular_mutation.tsv", sep = "\t", header = TRUE, data.table = FALSE)
df <- df[df$replicate == Sample,]

#Step 1: Define the lineage biased clones
if (Lineages[1] == "HSC"){
    stats <- read.table("../ChiSquare_statistics.tsv", sep = "\t", header = TRUE)
}else{
    stats <- read.table("../../PBMC/ChiSquare_statistics.tsv", sep = "\t", header = TRUE)
}
stats %>% arrange(pvalue) %>% subset(pvalue <= 0.05 & odds_ratio > 1) -> lin1
stats %>% arrange(pvalue) %>% subset(pvalue <= 0.05 & odds_ratio <= 1) -> lin2
lineages <- list("Lineage1" = lin1$mutation, "Lineage2" = lin2$mutation)
df <- annotate_clones(df, lineages)

#remove cells with multiple lineages
one <- df[df$bias == "Lineage1",]
two <- df[df$bias == "Lineage2",]
multi <- one[one$cell %in% two$cell,"cell"]
df <- df[!(df$cell %in% multi),]
print(table(df$bias))

#step 3: UMAP
SO <- annotate_SO(df, SO)
create_umap(SO, output_dir, Sample)

#step 4: Create density UMAPs
#https://github.com/TeamPerie/Cosgrove-et-al-2022/blob/main/Figure1/INTEGRATION/IntegrateDataset_Analysis.html
density_plot(SO, lineages = Lineages, output_dir, Sample)

#step 5:DEG between myeloid/lymphoid
genes <- deg_between_lineages(SO, output_dir, Sample=Sample)
fc_threshold = 0.1
p_threshold = 0.05

markers_up = genes[genes$avg_log2FC > fc_threshold & genes$p_val < p_threshold, ]
write(markers_up %>% rownames, paste0(output_dir,"/",Sample,"_up_gene.txt"))
markers_down = genes[genes$avg_log2FC < -(fc_threshold) & genes$p_val < p_threshold, ]
write(markers_down%>% rownames, paste0(output_dir,"/",Sample, "_down_gene.txt"))


