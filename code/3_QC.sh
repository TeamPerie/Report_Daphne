 #!/bin/bash
conda activate signac

samples='SRR11539031 SRR11539032 SRR11539033 SRR11539034'
for sample in $samples
do
	mkdir -p ../output/signac-analysis/${sample}
	cd ../output/signac-analysis/${sample}
	
	if [ ${sample} == "SRR11539031" ] || [ ${sample} == "SRR11539032" ] 
	then
		MinReads=25
	else
		MinReads=60
	fi
	
	Rscript ../../../code/Rscripts/QC.R $sample 1000 $MinReads 20
done




# MinPeaks=1000 #log10 >= 3
# MinReads=25 for BM and 60 for PBMC
# Depth=20




