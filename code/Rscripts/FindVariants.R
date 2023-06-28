library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly = TRUE)
samples <- strsplit(args[1], "_")[[1]]

SO1 <- read_rds(paste0("../",samples[1],"/SO_Annotated.rds"))
SO2 <- read_rds(paste0("../",samples[2],"/SO_Annotated.rds"))

SO <- merge(x=SO1, y=SO2, add.cell.ids = c(paste0(args[2],"1"), paste0(args[2],"2")), merge.data = TRUE)
write_rds(SO, "SO_merged.rds")

mito.data <-  read_rds(paste0("../",samples[1],"/mito.data.rds")) #reference alleles are equal between replicates
variable.sites <- IdentifyVariants(SO, assay = "mito", refallele = mito.data$refallele, low.coverage_threshold = 10)

png("variantplot.png")
p <- VariantPlot(variants = variable.sites, min.cells = 5, concordance.threshold = 0.5)
print(p)
dev.off()

#Filter mutations based on thresholds
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 & n_cells_conf_detected <= 1000 & strand_correlation >= 0.5 & vmr > 0.01
)
write.table(high.conf[,c(1,2,5)], "high.conf_table.txt", sep=",")

#Compute allele frequencies for each replicate separate
SO1 <- AlleleFreq(SO1, variants = high.conf$variant, assay = "mito")
DefaultAssay(SO1) <- "alleles"
write_rds(SO1, paste0("SO_", args[2],"1_variants.rds"))

SO2 <- AlleleFreq(SO2, variants = high.conf$variant, assay = "mito")
DefaultAssay(SO2) <- "alleles"
write_rds(SO1, paste0("SO_", args[2],"2_variants.rds"))


