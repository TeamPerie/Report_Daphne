#!/bin/bash

tissues='BM PBMC'
for tissue in $tissues
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${tissue}
	Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/MutantCells.R" ${tissue}
done
