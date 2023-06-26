print(getwd())
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(readr)
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
sample = args[1]

#Load Cellranger output
counts <- Read10X_h5(filename = paste0("../../../output/cellranger_atac_out/",sample,"/outs/filtered_peak_bc_matrix.h5"))
metadata <- read.csv(
    file = paste0("../../../output/cellranger_atac_out/",sample,"/outs/singlecell.csv"),
    header = TRUE,
    row.names = 1
)

#Load gene annotations from Ensembl and change to UCSC style
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

#Create Seurat object from chromatin assay of the cellranger output
#Only include peaks included in at least 10 cells
#Reference genome: Human GRCh38
assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    annotation = annotations,
    min.cells = 10,
    fragments = paste0("../../../output/cellranger_atac_out/",sample,"/outs/fragments.tsv.gz")
)
    
SO <- CreateSeuratObject(
    counts = assay,
    assay = 'peaks',
    genome = "GRCh38",
    meta.data = metadata
)

#QC metrics that were computed by cellranger-atac
SO$pct_reads_in_peaks <- SO$peak_region_fragments / SO$passed_filters * 100

#Compute TSS enrichment score and nucleosome banding pattern
SO <- TSSEnrichment(SO)
SO <- NucleosomeSignal(SO)

#Visualize QC metrics for each cell
png("QCplots.png")
p <- VlnPlot(SO, c("pct_reads_in_peaks",
             "passed_filters"), pt.size = 0.1, ncol = 2)
print(p)
dev.off()


png("QCplot.png")  
p <- ggplot(SO@meta.data, aes(x=pct_reads_in_peaks, y=passed_filters, 
                         color = passed_filters > as.numeric(args[2]) & 
                           pct_reads_in_peaks > as.numeric(args[3]))) +
  geom_point() +
  scale_y_continuous(trans='log10') +
  labs(x = "FRIP", y = "Unique Fragments") +
  theme_minimal() +
  geom_hline(yintercept = as.numeric(args[2]), color = "black", linetype = 2) +
  geom_vline(xintercept = as.numeric(args[3]), color = "black", linetype = 2) +
  scale_color_manual(
    values = c("black", "firebrick"),
    breaks = c(FALSE, TRUE)) +
  theme(legend.position = "none")
print(p)
dev.off()

#Filter out low quality cells
SO<- subset(
    x = SO,
    subset = passed_filters >= as.numeric(args[2]) &
      pct_reads_in_peaks >= as.numeric(args[3])
)

#Load mgatk output
mito.data <- ReadMGATK(dir = paste0("../../../output/mgatk_out/",sample,"/final/"))
write_rds(mito.data, "mito.data.rds")

#Subset to cells present in the scATAC-seq assat
mito <- CreateAssayObject(counts = mito.data$counts)
mito <- subset(mito, cells = colnames(SO))

#Add assay and metadata to the seurat object
SO[["mito"]] <- mito
SO <- AddMetaData(SO, metadata = mito.data$depth, col.name = "mtDNA_depth")

#Visualize mDNA sequencing depth
png("mtDepth_plot.png")
p <- VlnPlot(SO, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
print(p)
dev.off()

#Derived from Wenjie Sun
draw_depth_plot <- function(coverage, save_tsv = "", title = "") {
    names(coverage) = c("loci", "cell_id", "depth")
    d = coverage[, .(
        "depth 95%" = quantile(depth, 0.05) %>% as.numeric,
        "depth median" = quantile(depth, 0.50) %>% as.numeric,
        "depth 5%" = quantile(depth, 0.95) %>% as.numeric
        ), by = loci]
    d
    d_long = melt(d, id.vars="loci", measure.vars=c("depth 95%", "depth median", "depth 5%"), variable.name = "Quantile", value.name = "depth")
    g2 = ggplot(d_long) + aes(x = loci, y = depth, color = Quantile) + geom_line() + 
        coord_polar(theta = "x", start = 0, direction = 1, clip = "on") + 
        scale_y_log10() +
        labs(x = "chrM loc", y = "Depth", title = title) +
      theme(
        text = element_text(size = 20),
        plot.title = element_text(size = 24, face = "bold")
      )
    if (save_tsv != "") {
        write_tsv(d_long, save_tsv)
    }
    g2
}

d = fread(paste0("../../../output/mgatk_out/",sample,"/final/lib_id.coverage.txt.gz"))
d = d[d$V2 %in% colnames(SO), ]
png("Circular_mtDepth_before.png")
p1 <- draw_depth_plot(d, title = "Before Filtering")
print(p1)
dev.off()

#Filter out cells with low mtDNA depth
SO <- subset(SO, mtDNA_depth >= as.numeric(args[4]))
write_rds(SO, "SO_QCFiltered.rds")

png("Circular_mtDepth_after.png")
d = d[d$V2 %in% colnames(SO), ]
p2 <- draw_depth_plot(d, title = "After Filtering")
print(p2)
dev.off()





