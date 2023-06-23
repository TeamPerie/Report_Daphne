library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(readr)

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
SO$pct_reads_in_DNase <- SO$DNase_sensitive_region_fragments / SO$passed_filters * 100
SO$blacklist_ratio <- SO$blacklist_region_fragments / SO$peak_region_fragments

#Compute TSS enrichment score and nucleosome banding pattern
SO <- TSSEnrichment(SO)
SO <- NucleosomeSignal(SO)

#Visualize QC metrics for each cell
png("QCplots.png", width=1500, height=1000)
p <- VlnPlot(SO, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks", 
             "peak_region_fragments", "pct_reads_in_DNase", "blacklist_ratio"), pt.size = 0.1, ncol = 4)
print(p)
dev.off()

#Filter out low quality cells
SO<- subset(
    x = SO,
    subset = peak_region_fragments >= as.numeric(args[2]) &
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
png("depth_plot.png")
p <- VlnPlot(SO, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
print(p)
dev.off()

SO <- subset(SO, mtDNA_depth >= as.numeric(args[4]))
write_rds(SO, "SO_QCFiltered.rds")






