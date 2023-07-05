

volcano_plot <- function(genes, output_dir, Sample, fc_threshold, p_threshold, lineages, label){
    genes$diffexpressed <- "NO"
    genes$diffexpressed[genes$avg_log2FC > fc_threshold & genes$p_val < p_threshold] <- "UP"
    genes$diffexpressed[genes$avg_log2FC < -fc_threshold & genes$p_val < p_threshold] <- "DOWN"

    svg(paste0(output_dir, "/Volcano_", Sample, ".svg"))
    p <- ggplot(data = genes, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), col = "gray", linetype = 'dashed') +
        geom_hline(yintercept = -log10(p_threshold), col = "gray", linetype = 'dashed') + 
        geom_point()+
        scale_color_manual(values = c("NO" = "grey", 
                                    "UP" = ifelse(lineages[1] == "Lymphoid", "red", "green"),
                                    "DOWN" = ifelse(lineages[2] == "Myeloid", "blue", "black")),
                       labels = c(lineages[2], "Not significant", lineages[1]))  +
        geom_label(data = subset(genes, gene %in% label), 
          aes(x=avg_log2FC, y=-log10(p_val), label=gene), show_guide  = FALSE) +
        theme_classic()
    print(p)
    dev.off()
}


library(ggplot2)
library(dplyr)
library(data.table)
library(readr)

args = commandArgs(trailingOnly = TRUE)
Lineages <- strsplit(args, "_")[[1]]
output_dir = getwd()

genes1 <- read.table("Markers_BM1.tsv", sep = "\t", header = TRUE)
genes2 <- read.table("Markers_BM2.tsv", sep = "\t", header = TRUE)

up1 <- read_lines("BM1_up_gene.txt")
up2 <- read_lines("BM2_up_gene.txt")
down1 <- read_lines("BM1_down_gene.txt")
down2 <- read_lines("BM2_down_gene.txt")

up <- up1[up1 %in% up2]
down <- down1[down1 %in% down2]

fc_threshold = 0.1
p_threshold = 0.05

volcano_plot(genes1, output_dir, Sample="BM1", fc_threshold = fc_threshold, p_threshold = p_threshold, lineages = Lineages, 
    label = c(up,down))

volcano_plot(genes2, output_dir, Sample="BM2", fc_threshold = fc_threshold, p_threshold = p_threshold, lineages = Lineages, 
    label = c(up,down))