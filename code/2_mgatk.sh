#!/bin/bash

samples='SRR11539031 SRR11539032 SRR11539033 SRR11539034'
for sample in $samples
do
    mgatk tenx -i ../output/cellranger_atac_out/${sample}/outs/possorted_bam.bam \
      -g rCRS \
      -n ${sample} -o ../output/mgatk_out/${sample} -c 5 \
      -bt CB -b ../output/cellranger_atac_out/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
      --keep-duplicates \
      --alignment-quality 0
done




