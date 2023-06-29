#!/bin/bash

tissues='BM PBMC'
for tissue in $tissues
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	
	if [ ${tissue} == "BM"  ] 
	then
		lineages='HSC_Other'
	else
		lineages='Lymphoid_Myeloid'
	fi
	
	Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/LineageBias.R" ${lineages} ${tissue}
done
