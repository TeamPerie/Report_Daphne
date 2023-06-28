#!/bin/bash

tissues='BM PBMC'
for tissue in $tissues
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	
	if [ ${tissue} == "BM"  ] 
	then
		samples='SRR11539031_SRR11539032'
	else
		samples='SRR11539033_SRR11539034'
	fi
	
	Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/FindVariants.R" ${samples} ${tissue}
done
