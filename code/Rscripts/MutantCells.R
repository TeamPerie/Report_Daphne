library(Signac)
library(Seurat)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

get_mat = function(af_mat, bp_count_m) {
    mut_loc = rownames(af_mat)
    cell = colnames(af_mat)
    mut_loc_name = str_match(string = mut_loc, "(\\d+)([ATCG]).([ATCG])")
    fwd_name = paste0(mut_loc_name[, 4], "-", mut_loc_name[, 2], "-fwd")
    rev_name = paste0(mut_loc_name[, 4], "-", mut_loc_name[, 2], "-rev")

    fwd_mat = bp_count_m[fwd_name, cell]
    rev_mat = bp_count_m[rev_name, cell]

    list(af_mat = af_mat, fwd_mat = fwd_mat, rev_mat = rev_mat)
}

subset_mat = function(l, mut_v, cell_v) {
    i = rownames(l[[1]]) %in% mut_v
    lapply(l, function(x) {
        x[i, cell_v]
    })
}

preprocess_mut = function(SO, sample_name, af_threshold = 0.05, fwd_threshold = 2, rev_threshold = 2, max_cell_n = 1000) {

    d_meta = SO@meta.data

    af_mat <- as.data.frame(SO[['alleles']]@data)
    bp_count_m = SO[['mito']]@data

    ## Matrix list: af, fwd reads count, rev reads count
    mat_l = get_mat(af_mat, bp_count_m)

    mut_mat = data.matrix(mat_l$af_mat > af_threshold) & data.matrix(mat_l$fwd_mat >= fwd_threshold) & data.matrix(mat_l$rev_mat >= rev_threshold)

    n_cell = rowSums(mut_mat)
    print("Cell number per mutation:")
    table(n_cell) %>% print
    n_mut = colSums(mut_mat)
    print("Mutation number per cell:")
    table(n_mut) %>% print
    filtered_mut = names(n_cell)[n_cell < max_cell_n] 
    mat_l2 = subset_mat(mat_l, filtered_mut, rownames(d_meta))

    mut_mat2 = data.matrix(mat_l2$af_mat > af_threshold) & data.matrix(mat_l2$fwd_mat >= fwd_threshold) & data.matrix(mat_l2$rev_mat >= rev_threshold)
    n_cell2 = rowSums(mut_mat2)
    print("After filtering big clones:")
    print("Cell number per mutation:")
    table(n_cell2) %>% print
    n_mut2 = colSums(mut_mat2)
    print("Mutation number per cell:")
    table(n_mut2) %>% print


    mut_v = names(n_cell2)[n_cell2 > 0]

    d = lapply(mut_v, function(x) {
        x_i = which(rownames(mut_mat2) == x)
        y = which(mut_mat2[x_i, ])
        cell_barcode = colnames(mut_mat2)[y]
        data.table(
            cell = cell_barcode,
            mutation = x,
            af = mat_l2[[1]][x_i, y] %>% unlist,
            fwd = mat_l2[[2]][x_i, y] %>% unlist,
            rev = mat_l2[[3]][x_i, y] %>% unlist,
            annotation = d_meta[cell_barcode, "Cell.Types"],
            replicate = sample_name
        )
    }) %>% rbindlist
    d
}

repl = paste0(args[1],"1")
SO <- read_rds(paste0("SO_",repl,"_variants.rds"))
d1 = preprocess_mut(SO, repl, max_cell_n = 1000, af_threshold = 0.05, fwd_threshold = 2, rev_threshold = 2)

repl = paste0(args[1],"2")
SO <- read_rds(paste0("SO_",repl,"_variants.rds"))
d2 = preprocess_mut(SO, repl, max_cell_n = 1000, af_threshold = 0.05, fwd_threshold = 2, rev_threshold = 2)

d = rbind(d1, d2)

write_tsv(d, "table_cellular_mutation.tsv")

