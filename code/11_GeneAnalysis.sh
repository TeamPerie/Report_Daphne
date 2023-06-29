#!/bin/bash

lineages='HSC_Other Lymphoid_Myeloid'
#samples='SRR11539031 SRR11539032'
samples='BM1 BM2'
for lineage in $lineages
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/BM/${lineage}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/BM/${lineage}
	for sample in $samples
	do
		Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/DEG.R" ${lineage} ${sample}
	done
done
