
library(Signac)
library(BuenColors)
library(dplyr)
library(ggrastr)
library(SummarizedExperiment)

## Reproduce results Lareau
#https://github.com/caleblareau/mtscATACpaper_reproducibility
# Import processed data
cd34_clone_df <- readRDS("../../data/lareau/CD34_clone_DF.rds")
pbmc_clone_df <- readRDS("../../data/lareau/PBMC_clone_DF.rds")
cd34_mut_se <- readRDS("../../data/lareau/filteredCD34_mgatk_calls.rds")
pbmc_mut_se <- readRDS("../../data/lareau/filteredpbmcs_mgatk_calls.rds")

tb <- theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) 

make_4plot_grid <- function(variant){
  # Make plots of clones
  cd34_clone_df$color_AF <- assays(cd34_mut_se)[["allele_frequency"]][variant,]
  p_CD34_AF <- ggplot(cd34_clone_df %>% arrange(color_AF), aes(x= X1, y = X2, color = color_AF)) +
    geom_point_rast(size = 1) +
    tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) + ggtitle(paste0("BM ", variant))
  
  pbmc_clone_df$color_AF <- assays(pbmc_mut_se)[["allele_frequency"]][variant,]
  p_PBMC_AF <- ggplot(pbmc_clone_df %>% arrange(color_AF), aes(x= UMAP_1, y = UMAP_2, color = color_AF)) +
    geom_point_rast(size = 1) +
    tb + scale_color_gradientn(colors = c("lightgrey", "firebrick")) + ggtitle(paste0("PBMC ", variant))
  
  cowplot::plot_grid(p_CD34_AF,p_PBMC_AF, ncol = 2)
}

for (mut in c("12868G>A", "2788C>A", "3209A>G")){
    print(substr(mut, 1, -3))
	png(paste0("Lareau_VAF",substr(mut, 1, -3),".png"))
	p <- make_4plot_grid(mut)
	print(p)
	dev.off()
}

## Our results
plot_AF <- function(variant, pbmc1_af, pbmc1_umap, bm1_af, bm1_umap,
								pbmc2_af, pbmc2_umap, bm2_af, bm2_umap){
  PBMC1_AF <- data.frame(mut = pbmc1_af[,variant], 
                         rownames = rownames(pbmc1_af))
  df <- inner_join(PBMC1_AF, pbmc1_umap, by = "rownames") 
  df %>% arrange(mut) -> df
  p_PBMC1 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=as.numeric(mut))) +
    geom_point(size=1 +
    scale_color_gradient(low = "grey", high = "red")  +
    tb + ggtitle(paste0("PBMC1 ", variant))
  
  BM1_AF <- data.frame(mut = bm1_af[,variant], 
                         rownames = rownames(bm1_af))
  df <- inner_join(BM1_AF, bm1_umap, by = "rownames") 
  df %>% arrange(mut) -> df
  p_BM1 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=as.numeric(mut))) +
    geom_point(size=1 +
    scale_color_gradient(low = "grey", high = "red")  +
    tb + ggtitle(paste0("BM1 ", variant))
    
  PBMC2_AF <- data.frame(mut = pbmc2_af[,variant], 
                         rownames = rownames(pbmc2_af))
  df <- inner_join(PBMC2_AF, pbmc2_umap, by = "rownames") 
  df %>% arrange(mut) -> df
  p_PBMC2 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=as.numeric(mut))) +
    geom_point(size=1 +
    scale_color_gradient(low = "grey", high = "red")  +
    tb + ggtitle(paste0("PBMC2 ", variant))
  
  BM2_AF <- data.frame(mut = bm2_af[,variant], 
                         rownames = rownames(bm2_af))
  df <- inner_join(BM2_AF, bm2_umap, by = "rownames") 
  df %>% arrange(mut) -> df
  p_BM2 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=as.numeric(mut))) +
    geom_point(size=1 +
    scale_color_gradient(low = "grey", high = "red")  +
    tb + ggtitle(paste0("BM2 ", variant))
    
  cowplot::plot_grid(p_BM1, p_PBMC1, p_BM2, p_PBMC2, nrow = 2, ncol = 2)
}

PBMC1 <- readRDS("PBMC/SO_PBMC1_variants.rds")
PBMC1_AF <- as.data.frame(t(as.data.frame(as.matrix(PBMC1[["alleles"]]@data))))
PBMC1_UMAP <- as.data.frame(PBMC1@reductions$umap@cell.embeddings) 
PBMC1_UMAP %>% mutate(rownames = rownames(PBMC1_UMAP)) -> PBMC1_UMAP

BM1 <- readRDS("BM/SO_BM1_variants.rds")
BM1_AF <- as.data.frame(t(as.data.frame(as.matrix(BM1[["alleles"]]@data))))
BM1_UMAP <- as.data.frame(BM1@reductions$umap@cell.embeddings) 
BM1_UMAP %>% mutate(rownames = rownames(BM1_UMAP)) -> BM1_UMAP

PBMC2 <- readRDS("PBMC/SO_PBMC2_variants.rds")
PBMC2_AF <- as.data.frame(t(as.data.frame(as.matrix(PBMC2[["alleles"]]@data))))
PBMC2_UMAP <- as.data.frame(PBMC2@reductions$umap@cell.embeddings) 
PBMC2_UMAP %>% mutate(rownames = rownames(PBMC2_UMAP)) -> PBMC2_UMAP

BM2 <- readRDS("BM/SO_BM2_variants.rds")
BM2_AF <- as.data.frame(t(as.data.frame(as.matrix(BM2[["alleles"]]@data))))
BM2_UMAP <- as.data.frame(BM2@reductions$umap@cell.embeddings) 
BM2_UMAP %>% mutate(rownames = rownames(BM2_UMAP)) -> BM2_UMAP
  
for (mut in c("12868G>A", "2788C>A", "3209A>G")){
	png(paste0("Daphne_VAF",substr(mut, 1, -3),".png"))
	p <- plot_AF(mut,
						PBMC1_AF,
						PBMC1_UMAP,
						BM1_AF,
						BM1_UMAP,
						PBMC2_AF,
						PBMC2_UMAP,
						BM2_AF,
						BM2_UMAP)
	print(p)
	dev.off()
}