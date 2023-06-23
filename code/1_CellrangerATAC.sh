#!/bin/bash

cd ../output/cellranger_atac_out/

samples='SRR11539031 SRR11539032 SRR11539033 SRR11539034'
for sample in $samples
do
    cellranger-atac count \
      --id=${sample} \
      --fastqs=../data/fastq/${sample} \
      --reference=../data/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/ \
      --force-cells=6000 \
      --localcores=10 \
      --localmem=64
done