#!/bin/bash
#SBATCH --time=7:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=92768     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=CD_h2
#SBATCH --error=error/CD.error
#SBATCH --output=output/CD.output
#SBATCH --mail-user=brice7@unl.edu

module load plink

# Loop
for i in {1..1000}; do
  # Run PLINK command for each file
  plink --bfile combined_data --extract "Random_Marker_Sets_CD/RandomSet_${i}.txt" --make-bed --out "CobDiameter/data_${i}"

./ldak5.2.linux --calc-kins-direct "CobDiameter/kin_${i}" --bfile "CobDiameter/data_${i}" --power 1

./ldak5.2.linux --reml "CobDiameter/reml${i}" --pheno reseq.pheno.fam --mpheno 5 --grm "CobDiameter/kin_${i}"

rm CobDiameter/kin_${i}**
rm CobDiameter/data_${i}**
done