 #Torque Configuration
 #PBS -l walltime=24:00:00
 #PBS -l mem=24gb
 #PBS -l nodes=1:ppn=4
 #PBS -q batch
 #PBS -N Repro
 #PBS -e /data/users/mvanginn/CellLineageTracing/Reproduce_Thesis_InstitutCurie/log/DEG.err
 #PBS -o /data/users/mvanginn/CellLineageTracing/Reproduce_Thesis_InstitutCurie/log/DEG.log
 
 source ~/.bashrc
conda activate signac

cd "/data/users/mvanginn/CellLineageTracing/Reproduce_Thesis_InstitutCurie/code"
#bash /data/users/mvanginn/CellLineageTracing/Reproduce_Thesis_InstitutCurie/code/9_VAFplots.sh "/data/users/mvanginn/CellLineageTracing"

bash /data/users/mvanginn/CellLineageTracing/Reproduce_Thesis_InstitutCurie/code/11_GeneAnalysis.sh "/data/users/mvanginn/CellLineageTracing"