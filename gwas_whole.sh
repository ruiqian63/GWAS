#!/bin/bash

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=48:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=gwas_whole # Job name
#SBATCH --output=gwas_whole.%j.out # Name of the output file
#SBATCH --error=gwas_whole.%j.error
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of threads/tasks on one node
#SBATCH --mem=64G

## commands
source /programs/biogrids.shrc

# add phenotype CHD
awk '{print $2, 2}' case.fam > case_ids.txt
awk '{print $2, 1}' control.fam > control_ids.txt
cat case_ids.txt control_ids.txt > all_ids.txt
awk 'NR==FNR{a[$1]=$2; next} {if ($2 in a) $6=a[$2]; print}' all_ids.txt merged_whole.fam > updated_phenotype.fam

mv merged_whole.bim updated_phenotype.bim
mv merged_whole.bed updated_phenotype.bed

# quality control
plink --bfile  updated_phenotype --genome --min 0.1 --allow-no-sex --make-bed --out genome_filtered
plink --bfile genome_filtered --mind 0.1 --make-bed --out mind_filtered
plink --bfile mind_filtered --geno 0.02 --make-bed --out geno_filtered
plink --bfile geno_filtered --maf 0.05 --make-bed --out maf_filtered
plink --bfile maf_filtered --hwe 0.05 --make-bed --out qc_updated_phenotype_nosex_all

# pca
plink --bfile merged_data --pca 3 --out pca_results
awk '{print $1, $2, $3, $4, $5}' pca_results.eigenvec > pca_covar.txt

# logistic
plink --bfile qc_updated_phenotype_nosex_all --logistic --allow-no-sex --out gwas_logistic_results_whole
plink --bfile qc_updated_phenotype_nosex_all --logistic --covar pca_covar.txt --allow-no-sex --out gwas_logistic_results_whole_pca
