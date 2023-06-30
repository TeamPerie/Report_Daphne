lineages='Lymphoid_Myeloid HSC_Other'
samples='BM1 BM2'
for lineage in $lineages
do
	mkdir -p ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/BM/${lineage}
	cd ${1}/Reproduce_Thesis_InstitutCurie/output/signac_analysis/BM/${lineage}
	for sample in $samples
	do
		Rscript "${1}/Reproduce_Thesis_InstitutCurie/code/Rscripts/Pathways.R" ${lineage} ${sample}
	done
done
