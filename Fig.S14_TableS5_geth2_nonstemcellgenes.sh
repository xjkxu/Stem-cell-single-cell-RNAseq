#!/bin/bash
#SBATCH --time=5:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=92768     # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=h2
#SBATCH --error=GeneSet6/h2.error
#SBATCH --output=GeneSet6/h2.output
#SBATCH --mail-user=brice7@unl.edu

module load plink 

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/CD_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_CD --pheno reseq.pheno.fam --mpheno 1 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/CW_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_CW --pheno reseq.pheno.fam --mpheno 2 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/ED_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_ED --pheno reseq.pheno.fam --mpheno 3 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/EL_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_EL --pheno reseq.pheno.fam --mpheno 4 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/SSL_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_SSL --pheno reseq.pheno.fam --mpheno 5 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/ERN_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_ERN --pheno reseq.pheno.fam --mpheno 6 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/EW_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_EW --pheno reseq.pheno.fam --mpheno 7 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/ERkN_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_ERkN --pheno reseq.pheno.fam --mpheno 8 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**

#Get h2 for these marker setsplink --bfile fullgenome --extract GeneSet6/ColocalizedMarkers/KW20_ColocalizedMarkers.txt --make-bed --out GeneSet6/ColocalizedMarkers/data./ldak5.2.linux --calc-kins-direct GeneSet6/ColocalizedMarkers/kin --bfile GeneSet6/ColocalizedMarkers/data --power 1./ldak5.2.linux --reml GeneSet6/ColocalizedMarkers/reml_KW20 --pheno reseq.pheno.fam --mpheno 9 --grm GeneSet6/ColocalizedMarkers/kinrm GeneSet6/ColocalizedMarkers/kin**rm GeneSet6/ColocalizedMarkers/data**