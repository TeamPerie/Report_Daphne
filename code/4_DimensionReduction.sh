 #!/bin/bash

samples='SRR11539031 SRR11539032 SRR11539033 SRR11539034'
for sample in $samples
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${sample}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/${sample}
	
	if [ ${sample} == "SRR11539031" ] || [ ${sample} == "SRR11539032" ] 
	then
		tissue="BM"
	else
		tissue="PBMC"
	fi
	
	Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/DimensionReduction.R" $tissue
done
