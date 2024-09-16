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

# convert .vcf files
plink --vcf merged_final_cleaned_no_multiple_rsids.vcf.gz --make-bed --out merged_whole_nomulti

# add phenotype CHD
awk '{print $2, 2}' case.fam > case_ids.txt
awk '{print $2, 1}' control.fam > control_ids.txt
cat case_ids.txt control_ids.txt > all_ids.txt
awk 'NR==FNR{a[$1]=$2; next} {if ($2 in a) $6=a[$2]; print}' all_ids.txt merged_whole_nomulti.fam > updated_phenotype_nomulti.fam

mv merged_whole_nomulti.bim updated_phenotype_nomulti.bim
mv merged_whole_nomulti.bed updated_phenotype_nomulti.bed

# quality control
plink --bfile  updated_phenotype_nomulti --genome --min 0.05 --allow-no-sex --make-bed --out genome_filtered_nomulti
plink --bfile genome_filtered_nomulti --mind 0.05 --make-bed --out mind_filtered_nomulti
plink --bfile mind_filtered_nomulti --geno 0.02 --make-bed --out geno_filtered_nomulti
plink --bfile geno_filtered_nomulti --maf 0.05 --make-bed --out maf_filtered_nomulti
plink --bfile maf_filtered_nomulti --hwe 0.05 --make-bed --out qc_updated_phenotype_nosex_all_nomulti

# pca
plink --bfile qc_updated_phenotype_nosex_all_nomulti --pca 3 --out pca_results_nomulti
awk '{print $1, $2, $3, $4, $5}' pca_results_nomulti.eigenvec > pca_covar_nomulti.txt

# logistic
plink --bfile qc_updated_phenotype_nosex_all_nomulti --logistic --allow-no-sex --out gwas_logistic_results_whole_nomulti
plink --bfile qc_updated_phenotype_nosex_all_nomulti --logistic --covar pca_covar_nomulti.txt --allow-no-sex --out gwas_logistic_results_whole_pca_nomulti

