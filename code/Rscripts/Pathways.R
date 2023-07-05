#Find enriched pathways
# Go to https://david.ncifcrf.gov/tools.jsp
# Step 1: upload the up/down regulated genes from BM*_up/down_genes.txt
# Step 2: OFFICIAL_GENE_SYMBOL & Homo sapiens
# Step 3: Gene List
# Go to KEGG_PATHWAY and download the file as BM*_up/down_KEGG.txt


pathway_plot <- function(df, lineage, output_dir, Sample){
    if(lineage == "Myeloid"){color = "blue"}
    if(lineage == "Lymphoid"){color = "red"}
    if(lineage == "HSC"){color = "green"}
    if(lineage == "Other"){color = "black"}
    df %>% separate(Term, into = c("hsa", "Term"), sep = ":", remove = TRUE) -> df
    df$Term <- factor(df$Term, levels=df$Term[order(df$PValue)])
    svg(paste0(output_dir, "/",Sample,"_Pathway_",lineage,".svg"))
    p<-ggplot(data=df, aes(x=Term, y=-log10(PValue))) +
      geom_bar(stat="identity", fill = color) + coord_flip() + 
      labs(title = lineage) + theme_classic() +
      theme(text = element_text(size = 14),
    		axis.title = element_text(size = 14),
    		axis.text = element_text(size = 14))

    print(p)
    dev.off()
}

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly = TRUE)
Lineages <- strsplit(args, "_")[[1]]
output_dir = getwd()
Sample = args[2]
print(Sample)

p_up <- read.table(paste0(output_dir, "/", Sample, "_up_KEGG.txt"), sep = "\t", header = TRUE)
pathway_plot(p_up, Lineages[1], output_dir, Sample)

p_down <- read.table(paste0(output_dir, "/", Sample, "_down_KEGG.txt"), sep = "\t", header = TRUE)
pathway_plot(p_down, Lineages[2], output_dir, Sample)