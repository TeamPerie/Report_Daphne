library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

#pbmc
pbmc <- fread("PBMC/table_cellular_mutation.tsv", sep = "\t", header = TRUE, data.table = FALSE)
pbmc %>% group_by(mutation, replicate) %>% summarize(count=n()) %>% 
spread(key = replicate, value = count, fill=0)-> clone_size  
r <- cor.test(clone_size$PBMC1, clone_size$PBMC2)
svg("PBMC/CloneSize.svg")
p <- ggplot(clone_size, aes(x=PBMC1, y=PBMC2)) +
    geom_point(aes(color = ifelse(PBMC1 == 0 | PBMC2 == 0, "A", "B")), size=2, show.legend = FALSE) +
    scale_color_manual(values = c("blue", "firebrick")) +
    scale_x_continuous(trans="log10") + 
    scale_y_continuous(trans="log10") + 
    labs( title = "Clone size correlation PBMC") +
    annotate("text", x = 5, y = 120, label = paste0("R = ", round(unname(r$estimate), digits=2))) +
    annotate("text", x = 5, y = 100, label = paste0("p = ", round(r$p.value,, digits=2))) +
    theme_classic() + theme(legend.position = "none", text = element_text(size = 22),
    		axis.title = element_text(size = 22),
    		axis.text = element_text(size = 16))
print(p)
dev.off()

#bm
bm<- fread("BM/table_cellular_mutation.tsv", sep = "\t", header = TRUE, data.table = FALSE)
bm %>% group_by(mutation, replicate) %>% summarize(count=n()) %>% 
spread(key = replicate, value = count, fill=0)-> clone_size  
r <- cor.test(clone_size$BM1, clone_size$BM2)
svg("BM/CloneSize.svg")
p <- ggplot(clone_size, aes(x=BM1, y=BM2)) +
    geom_point(aes(color = ifelse(BM1 == 0 | BM2 == 0, "A", "B")), size=2, show.legend = FALSE) +
    scale_color_manual(values = c("blue", "firebrick")) +
    scale_x_continuous(trans="log10") + 
    scale_y_continuous(trans="log10") + 
    labs( title = "Clone size correlation BM") +
    annotate("text", x = 5, y = 120, label = paste0("R = ", round(unname(r$estimate), digits=2))) +
    annotate("text", x = 5, y = 100, label = paste0("p = ", round(r$p.value,, digits=2))) +
    theme_classic() + theme(legend.position = "none", text = element_text(size = 22),
    		axis.title = element_text(size = 22),
    		axis.text = element_text(size = 16))
print(p)
dev.off()

#BM vs PBMC
pbmc$Sample <- rep("PBMC", nrow(pbmc))
bm$Sample <- rep("BM", nrow(bm))
df <- rbind(bm, pbmc)
df %>% group_by(mutation, Sample) %>% summarize(count=n()) %>% 
spread(key = Sample, value = count, fill=0)-> clone_size  
r <- cor.test(clone_size$BM, clone_size$PBMC)
svg("CloneSize.svg")
p <- ggplot(clone_size, aes(x=BM, y=PBMC)) +
    geom_point(aes(color = ifelse(BM == 0 | PBMC == 0, "A", "B")), size=2, show.legend = FALSE) +
     scale_color_manual(values = c("blue", "firebrick")) +
    scale_x_continuous(trans="log10") + 
    scale_y_continuous(trans="log10") + 
    labs( title = "Clone size correlation BM/PBMC") +
    annotate("text", x = 5, y = 120, label = paste0("R = ", round(unname(r$estimate), digits=2))) +
    annotate("text", x = 5, y = 100, label = "p = 1.61e-33") +
    theme_classic() + theme(legend.position = "none", text = element_text(size = 22),
    		axis.title = element_text(size = 22),
    		axis.text = element_text(size = 16))
print(p)
dev.off()


