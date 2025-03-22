#!/bin/bash
#SBATCH --time=5:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=92768     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=h2
#SBATCH --error=h2.error
#SBATCH --output=h2.output
#SBATCH --mail-user=brice7@unl.edu

module load plink

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/CD_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_CD" --pheno ../reseq.pheno.fam --mpheno 1 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/CW_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_CW" --pheno ../reseq.pheno.fam --mpheno 2 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/ED_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_ED" --pheno ../reseq.pheno.fam --mpheno 3 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/EL_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_EL" --pheno ../reseq.pheno.fam --mpheno 4 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/SSL_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_SSL" --pheno ../reseq.pheno.fam --mpheno 5 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/ERN_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_ERN" --pheno ../reseq.pheno.fam --mpheno 6 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/EW_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_EW" --pheno ../reseq.pheno.fam --mpheno 7 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/ERkN_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_ERkN" --pheno ../reseq.pheno.fam --mpheno 8 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile ../fullgenome --extract "ColocalizedMarkers/KW20_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"./../ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --allow-multi YES --power 1./../ldak5.2.linux --reml "ColocalizedMarkers/reml_KW20" --pheno ../reseq.pheno.fam --mpheno 9 --grm "ColocalizedMarkers/kin"rm ColocalizedMarkers/kin**rm ColocalizedMarkers/data**