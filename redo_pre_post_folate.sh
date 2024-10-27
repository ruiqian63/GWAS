#!/bin/bash

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=48:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=folate_gwas # Job name
#SBATCH --output=folate_gwas.%j.out # Name of the output file
#SBATCH --error=folate_gwas.%j.error
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of threads/tasks on one node
#SBATCH --mem=64G

## commands
source /programs/biogrids.shrc

#awk '$3 == "TRUE"' folate_cohort.txt > post_folate.txt
#awk '$3 == "FALSE"' folate_cohort.txt > pre_folate.txt
#cut -f1 post_folate.txt > post_folate_iids.txt
#cut -f1 pre_folate.txt > pre_folate_iids.txt

#filter the merged cohort by folate
bcftools view -S post_folate_iids.txt -Oz -o post_folate.vcf.gz merged_folate.vcf.gz
bcftools view -S pre_folate_iids.txt -Oz -o pre_folate.vcf.gz merged_folate.vcf.gz

bcftools index pre_folate.vcf.gz
bcftools index post_folate.vcf.gz

plink --vcf pre_folate.vcf.gz --make-bed --out pre_folate
plink --vcf post_folate.vcf.gz --make-bed --out post_folate

#construct phenotype_list
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' pre_folate.txt > pre_phenotype_map.txt
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' post_folate.txt > post_phenotype_map.txt

#add phenotype
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' pre_phenotype_map.txt pre_folate.fam > pre_updated_phenotype.fam
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' post_phenotype_map.txt post_folate.fam > post_updated_phenotype.fam
mv pre_folate.bim pre_updated_phenotype.bim
mv pre_folate.bed pre_updated_phenotype.bed
mv post_folate.bim post_updated_phenotype.bim
mv post_folate.bed post_updated_phenotype.bed

plink --bfile pre_updated_phenotype --mind 0.1 --make-bed --out pre_mind_filtered
plink --bfile pre_mind_filtered --geno 0.1 --make-bed --out pre_geno_filtered
plink --bfile pre_geno_filtered --maf 0.05 --make-bed --out pre_maf_filtered
plink --bfile pre_maf_filtered --hwe 0.05 --make-bed --out pre_qc_updated_phenotype_all

plink --bfile post_updated_phenotype --mind 0.1 --make-bed --out post_mind_filtered
plink --bfile post_mind_filtered --geno 0.1 --make-bed --out post_geno_filtered
plink --bfile post_geno_filtered --maf 0.05 --make-bed --out post_maf_filtered
plink --bfile post_maf_filtered --hwe 0.05 --make-bed --out post_qc_updated_phenotype_all

#plink --bfile pre_qc_updated_phenotype_all --pca 10 --out pre_pca_results
#awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' pre_pca_results.eigenvec > pre_pca_covar_10.txt
#plink --bfile post_qc_updated_phenotype_all --pca 10 --out post_pca_results
#awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' post_pca_results.eigenvec > post_pca_covar_10.txt

plink --bfile pre_qc_updated_phenotype_all --pca 3 --out pre_pca_results
awk '{print $1, $2, $3, $4, $5}' pre_pca_results.eigenvec > pre_pca_covar_3.txt
plink --bfile post_qc_updated_phenotype_all --pca 3 --out post_pca_results
awk '{print $1, $2, $3, $4, $5}' post_pca_results.eigenvec > post_pca_covar_3.txt



#plink2 --bfile pre_qc_updated_phenotype_all --glm --covar pre_pca_covar_10.txt --covar-variance-standardize --allow-no-sex --out pre_gwas_logistic_results_all_pca10_plink2 
#plink2 --bfile pre_qc_updated_phenotype_all --glm --covar pre_pca_covar_10.txt --covar-variance-standardize --allow-no-sex --out pre_gwas_logistic_results_all_pca10_plink2 --adjust
#plink2 --bfile post_qc_updated_phenotype_all --glm --covar post_pca_covar_10.txt --covar-variance-standardize --allow-no-sex --out post_gwas_logistic_results_all_pca10_plink2
#plink2 --bfile post_qc_updated_phenotype_all --glm --covar post_pca_covar_10.txt --covar-variance-standardize --allow-no-sex --out post_gwas_logistic_results_all_pca10_plink2 --adjust
