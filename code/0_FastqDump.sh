#!/bin/bash

samples='SRR11539031 SRR11539032 SRR11539033 SRR11539034'
for sample in $samples
do
	mkdir -p ../data/fastq/${sample}
	cd ../data/fastq/${sample}
	fastq-dump --origfmt -I --split-files --gzip ../data/sra/${sample}/${sample}.sra
done
