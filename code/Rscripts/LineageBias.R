
##Lineage bias analysis with the chisquare test
#
# df = dataframe containing cells, mutations and annotations
# lineage = the lineages to test
# output_dir = directory to save results
##
chisquare <- function(
    df,
    lineages,
    output_dir,
    Sample,
    remove_other = FALSE
){
  df <- annotate_lineages(df, lineages, remove_other)
  create_umap(df, output_dir, Sample, lineages)
  clones <- stacked_barplot(df, output_dir, lineages=lineages)
  df2 <- df[df$mutation %in% clones,]
  stats <- chisquare_stats(df2, output_dir, lineages)
  ranked_clones(stats, output_dir, lineages)
  df_sign <- stats %>% subset(pvalue < 0.05)
  heatmap_lineages(df_sign, output_dir, lineages)
}

create_umap <- function(df, output_dir, Sample, lineages){
  for (repl in paste0(Sample, c(1,2))){
    SO <- readRDS(paste0("SO_",repl,"_variants.rds"))
    SO_UMAP <- as.data.frame(SO@reductions$umap@cell.embeddings)
    SO_UMAP %>% mutate(cell = rownames(SO_UMAP)) -> SO_UMAP
    df2 <- inner_join(SO_UMAP, df[,c("cell","annotation")], by = "cell") 
    if (lineages[1]=="HSC"){colors <- c("HSC"="green","Other"="black")}else{
      colors <- c("Lymphoid" = "red", "Myeloid" = "blue", "NK-Cell" = "grey", "DC" = "grey")}
    svg(paste0(output_dir,"/", repl, "_UMAP_lineages.svg"))
    p <- ggplot(df2, aes(x = UMAP_1, y = UMAP_2, color = annotation)) +
      geom_point() +
      scale_color_manual(values = colors) +
      theme_classic() + theme(legend.position = "none")
    print(p)
    dev.off()
  }
}


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

## Create a stacked barplot with the proportion of cells per clone in each lineage
#
# df = dataframe containing cells, mutations and annotations
# output_dir = directory to save results
# lineages = between which lineages?
##
stacked_barplot <- function(
    df,
    output_dir,
    lineages = NULL
){
  df %>% group_by(mutation, annotation) %>% summarize(count=n()) %>%
    ungroup() %>%
    group_by(mutation) %>%
    mutate(total = sum(count)) %>%
    mutate(proportion = count/total) %>%
    dplyr::filter(total >= 5) %>% subset(annotation %in% lineages)  %>% 
    group_by(annotation) %>%
    arrange(annotation, desc(proportion)) -> df2
  df2$mutation <- factor(df2$mutation, levels = unique(df2$mutation))
  df2$all <- "all"
  
  svg(paste0(output_dir, "/Stacked_barplot_clones.svg"))
  p1 <- ggplot(df2, aes(x=mutation, fill=annotation, y=proportion)) + 
    geom_bar(width=1, position = "fill", stat = "identity") + theme_classic() +
    scale_fill_manual(values = c("HSC" = "green", "Lymphoid" = "red","Other"= "black", "Myeloid" = "blue")) +
    theme(text = element_text(size = 14), axis.title = element_text(size = 14), axis.text = element_text(size = 14),
	axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
	axis.ticks.y = element_blank()) + 
    labs(x = paste0("Clones (n=", length(unique(df2$mutation)), ")"), y = element_blank())
  p2 <- ggplot(df2, aes(x=all, fill=annotation, y=proportion)) + 
    scale_fill_manual(values = c("HSC" = "green", "Lymphoid" = "red","Other"= "black", "Myeloid" = "blue")) +
    geom_bar(width=1, position = "fill", stat = "identity") + theme_classic() +
    theme(text = element_text(size = 14), axis.title = element_text(size = 14),
	axis.text = element_text(size = 14), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
	legend.position = "none")
  print(grid.arrange(p2, p1, ncol = 2, widths = c(1,3)))
  dev.off()
  
  return(levels(df2$mutation))
}

## Perform the ChiSquare test
#
# df = dataframe containing cells, mutations and annotations
# output_dir = directory to save results
##
chisquare_stats <- function(
    df,
    output_dir,
    lineages
){
  #Vector with the total number of cells per lineage in the entire dataset
  df <- df[df$annotation %in% lineages,]
  total_df <- df %>% group_by(annotation) %>% 
    summarize(count = n()) %>% 
    mutate(total = sum(count)) %>%
    mutate(proportion = count/total)
  total_vector <- total_df$proportion
  names(total_vector) <- total_df$annotation
  
  #Test for each clone
  if (all(sort(unique(df$annotation)) == sort(c("Lymphoid", "Myeloid")))){
    cs_df <- data.frame(mutation = NA, stat = NA, pvalue = NA, n_cells = NA, perc_Lymphoid = NA, perc_Myeloid= NA, odds_ratio = NA)
    for (x in unique(df$mutation)){
      clone_df <- df %>% filter(mutation == x) %>% group_by(annotation) %>% summarize(count = n())
      clone_vector <- clone_df$count
      names(clone_vector) <- clone_df$annotation
      clone_vector <- unname(clone_vector[sort(unique(total_df$annotation))])
      clone_vector <- ifelse(is.na(clone_vector), 0, clone_vector)
      stat <- chisq.test(clone_vector, p = total_vector)
      
      #normalize for clone size
#       perc_Lymphoid = clone_vector[1]/sum(clone_vector)
#       perc_Myeloid = clone_vector[2]/sum(clone_vector)
    
      #normalize for lineage size
      perc_Lymphoid = clone_vector[1]/total_df$count[1]
      perc_Myeloid = clone_vector[2]/total_df$count[2]
      stat_df <- data.frame(mutation = x, stat = unname(stat$statistic), pvalue = unname(stat$p.value), 
                            n_cells = sum(clone_df$count),
                            perc_Lymphoid = perc_Lymphoid,
                            perc_Myeloid = perc_Myeloid,
                            odds_ratio = perc_Lymphoid/perc_Myeloid)
      cs_df <- rbind(cs_df, stat_df)
    }
  }
  if (all(sort(unique(df$annotation)) == sort(c("HSC", "Other")))){
    cs_df <- data.frame(mutation = NA, stat = NA, pvalue = NA, n_cells = NA, perc_HSC= NA, perc_Other= NA, odds_ratio = NA)
    for (x in unique(df$mutation)){
      clone_df <- df %>% filter(mutation == x) %>% group_by(annotation) %>% summarize(count = n())
      clone_vector <- clone_df$count
      names(clone_vector) <- clone_df$annotation
      clone_vector <- unname(clone_vector[sort(unique(total_df$annotation))])
      clone_vector <- ifelse(is.na(clone_vector), 0, clone_vector)
      stat <- chisq.test(clone_vector, p = total_vector)
      
      #normalize for clone size
#       perc_Lymphoid = clone_vector[1]/sum(clone_vector)
#       perc_Myeloid = clone_vector[2]/sum(clone_vector)
    
      #normalize for lineage size
      perc_HSC = clone_vector[1]/total_df$count[1]
      perc_Other = clone_vector[2]/total_df$count[2]
      stat_df <- data.frame(mutation = x, stat = unname(stat$statistic), pvalue = unname(stat$p.value), 
                            n_cells = sum(clone_df$count),
                            perc_HSC = perc_HSC,
                            perc_Other = perc_Other,
                            odds_ratio = perc_HSC/perc_Other)
      cs_df <- rbind(cs_df, stat_df)
    }
  }  
  cs_df %<>% drop_na() %>% arrange(pvalue)
  cs_df$p_adjust <- p.adjust(cs_df$pvalue, "bonferroni")
  write.table(cs_df, file=paste0(output_dir, "/ChiSquare_statistics.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  return(cs_df)
}

## Create a plot of the lineage biased ranked clones
#
# df = dataframe containing mutations lineage bias p-values
# output_dir = directory to save results
##
ranked_clones <- function(df, output_dir, lineages){
  df %>% arrange(pvalue) -> stats
  stats$rank <- rownames(stats)
  stats$Bias <- "unbiased"
  stats[stats$odds_ratio > 1 & stats$pvalue < 0.05,"Bias"] <- lineages[1]
  stats[stats$odds_ratio < 1 & stats$pvalue < 0.05,"Bias"] <- lineages[2]
  
  svg(paste0(output_dir, "/RankedClones.svg"))
  p <- ggplot(stats, aes(x=as.numeric(rank), y=-log10(pvalue), colour = as.factor(Bias))) +
    geom_point() +
    scale_color_manual("Biased clones", values = c(ifelse(lineages[1] == "Lymphoid", "red", "green"),
	ifelse(lineages[2] == "Myeloid", "blue", "black"), "grey")) +
    geom_text(data = subset(stats, pvalue < 0.05), 
              aes(x=as.numeric(rank), y=-log10(pvalue), label=mutation),
              hjust = 0, nudge_x = 3) +
    labs(x = "Rank sorted clones", y = "-log10 p-value") + theme_classic() +
    theme(text = element_text(size = 14),
    		axis.title = element_text(size = 14),
    		axis.text = element_text(size = 14))

  print(p)
  dev.off()
}

## Create a heatmap of the percentage of cells in each lineage per clone
#
# df = dataframe containing mutations and the proportion of cells per lineage
# output_dir = directory to save results
##
heatmap_lineages <- function(df, output_dir, lineages){
  if (all(sort(lineages) == sort(c("Lymphoid", "Myeloid")))){
#     df$perc_Lymphoid <- lapply(df_sign$perc_Lymphoid, asinh)
#     df$perc_Myeloid <- lapply(df_sign$perc_Myeloid, asinh)

    df$Lymphoid <- log10(df$perc_Lymphoid + 0.0001)
    df$Myeloid <- log10(df$perc_Myeloid + 0.0001)
    long_df <- gather(df, lineage, proportion, Lymphoid:Myeloid, factor_key=TRUE)
    dendro <- as.dendrogram(hclust(d=dist(x=df[,c("Lymphoid", "Myeloid")])))
    df_order <- order.dendrogram(dendro)
    long_df$mutation <- factor(x=long_df$mutation, levels = df$mutation[df_order], ordered=TRUE)
  }
  if (all(sort(lineages) == sort(c("HSC", "Other")))){
    df$HSC<- log10(df$perc_HSC + 0.0001)
    df$Other <- log10(df$perc_Other + 0.0001)
    long_df <- gather(df, lineage, proportion, HSC:Other, factor_key=TRUE)
    dendro <- as.dendrogram(hclust(d=dist(x=df[,c("HSC", "Other")])))
    df_order <- order.dendrogram(dendro)
    long_df$mutation <- factor(x=long_df$mutation, levels = df$mutation[df_order], ordered=TRUE)
  }
  
  svg(paste0(output_dir, "/Heatmap_LineageBias.svg"))
  plot <- ggplot(long_df, aes(y = mutation, x = lineage, fill = as.numeric(proportion))) + 
    geom_tile() + pretty_plot(fontsize = 15) + L_border() +
    scale_fill_gradientn("Proportion mutant cells", colors = jdb_palette("solar_rojos")) +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
  print(plot)
  dev.off()
}



library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(BuenColors)
library(stringr)
library(gridExtra)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
lineages <- strsplit(args[1], "_")[[1]]
df <- fread("table_cellular_mutation.tsv", sep = "\t", header = TRUE, data.table = FALSE)
#df$annotation <- as.character(df$annotation)

chisquare(df = df, lineages = lineages, remove_other = TRUE, output_dir = getwd(), Sample = args[2])

