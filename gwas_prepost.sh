#!/bin/bash

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=48:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=pre&post # Job name
#SBATCH --output=pre&post.%j.out # Name of the output file
#SBATCH --error=pre&post.%j.error
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of threads/tasks on one node
#SBATCH --mem=48G

## commands
source /programs/biogrids.shrc

#get vcf index
bcftools index post_folate_whole.vcf.gz 
bcftools index pre_folate_whole.vcf.gz

# convert .vcf files
plink --vcf post_folate_whole.vcf.gz --make-bed --out post_whole
plink --vcf pre_folate_whole.vcf.gz --make-bed --out pre_whole

#construct phenotype_list
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' pre_folate.txt > pre_phenotype_map.txt
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' post_folate.txt > post_phenotype_map.txt

#add phenotype
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' pre_phenotype_map.txt pre_whole.fam > pre_updated_phenotype.fam
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' post_phenotype_map.txt post_whole.fam > post_updated_phenotype.fam

mv pre_whole.bim pre_updated_phenotype.bim
mv pre_whole.bed pre_updated_phenotype.bed
mv post_whole.bim post_updated_phenotype.bim
mv post_whole.bed post_updated_phenotype.bed

#quality contory
plink --bfile pre_updated_phenotype --genome --min 0.05 --allow-no-sex --make-bed --out pre_genome_filtered
plink --bfile pre_genome_filtered --mind 0.05 --make-bed --out pre_mind_filtered
plink --bfile pre_mind_filtered --geno 0.02 --make-bed --out pre_geno_filtered
plink --bfile pre_geno_filtered --maf 0.05 --make-bed --out pre_maf_filtered
plink --bfile pre_maf_filtered --hwe 0.05 --make-bed --out pre_qc_updated_phenotype

#quality contory
plink --bfile post_updated_phenotype --genome --min 0.05 --allow-no-sex --make-bed --out post_genome_filtered
plink --bfile post_genome_filtered --mind 0.05 --make-bed --out post_mind_filtered
plink --bfile post_mind_filtered --geno 0.02 --make-bed --out post_geno_filtered
plink --bfile post_geno_filtered --maf 0.05 --make-bed --out post_maf_filtered
plink --bfile post_maf_filtered --hwe 0.05 --make-bed --out post_qc_updated_phenotype

#pca
plink --bfile pre_qc_updated_phenotype --pca 10 --out pre_pca_results
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' pre_pca_results.eigenvec > pre_pca_covar.txt
plink --bfile post_qc_updated_phenotype --pca 10 --out post_pca_results
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' post_pca_results.eigenvec > post_pca_covar.txt


#GWAS
plink --bfile pre_qc_updated_phenotype --logistic --allow-no-sex --out gwas_logistic_results_pre_folate
plink --bfile post_qc_updated_phenotype --logistic --allow-no-sex --out gwas_logistic_results_post_folate
plink --bfile pre_qc_updated_phenotype --logistic --covar pre_pca_covar.txt --allow-no-sex --out gwas_logistic_results_pre_folate_pca
plink --bfile post_qc_updated_phenotype --logistic --covar post_pca_covar.txt --allow-no-sex --out gwas_logistic_results_post_folate_pca


