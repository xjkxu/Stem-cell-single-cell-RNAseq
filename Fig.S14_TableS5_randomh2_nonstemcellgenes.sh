#!/bin/bash
#SBATCH --time=7:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=92768     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=ED_h2
#SBATCH --error=error/ED.error
#SBATCH --output=output/ED.output
#SBATCH --mail-user=brice7@unl.edu

module load plink
mkdir ED

# Loop
for i in {1..1000}; do
  # Run PLINK command for each file
  plink --bfile ../../fullgenome --extract "Random_Marker_Sets_ED/RandomSet_${i}.txt" --make-bed --out "ED/data_${i}"

./../../ldak5.2.linux --calc-kins-direct "ED/kin_${i}" --bfile "ED/data_${i}" --power 1

./../../ldak5.2.linux --reml "ED/reml${i}" --pheno ../../reseq.pheno.fam --mpheno 8 --grm "ED/kin_${i}"

rm ED/kin_${i}**
rm ED/data_${i}**
done