#!/bin/bash

#compress vcf
bgzip file.vcf

#unzip vcf.gz
gunzip file.vcf.gz

#delete multiple rsIDs
bcftools view --types snps -Oz -o out.vcf.gz input.vcf.gz

#get vcf index
bcftools index qc_control.vcf.gz
bcftools index qc_case.vcf.gz

# finding shared variants
bcftools query -f '%CHROM\t%POS\n' folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz | sort -u > file1_positions.txt
bcftools query -f '%CHROM\t%POS\n' folate_pass_pcgc_ctd_srv.vcf.gz | sort -u > file2_positions.txt
comm -12 <(sort file1_positions.txt) <(sort file2_positions.txt) > shared_positions.txt
vcftools --gzvcf folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz --positions shared_positions.txt --recode --out filtered_folate_pass_coding_sfari_sibs
vcftools --gzvcf folate_pass_pcgc_ctd_srv.vcf.gz --positions shared_positions.txt --recode --out filtered_folate_pass_pcgc_ctd_srv

# change vcf header
bcftools view -h folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz > hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz > folate_pass_coding_sfari_sibs_removed_RG.vcf.gz

# merge vcf
bcftools merge qc_control.vcf.gz qc_case.vcf.gz -O z -o merged_filtered.vcf.gz

# convert .vcf files
plink --vcf input_file.vcf.gz --make-bed --out output_prefix

# add sex information
awk '{if(NR>1) print $1, $2, $4}' file.sexcheck > updated_sex.txt
plink --bfile original_data --update-sex updated_sex.txt --make-bed --out updated_data

# add phenotype CHD
awk '{print $2, 2}' case.fam > case_ids.txt
awk '{print $2, 1}' control.fam > control_ids.txt
cat case_ids.txt control_ids.txt > all_ids.txt
awk 'NR==FNR{a[$1]=$2; next} {if ($2 in a) $6=a[$2]; print}' all_ids.txt updated_sex_all.fam > updated_phenotype.fam

# quality control
plink --bfile all_sex_phenotype --genome --min 0.1 --allow-no-sex --make-bed --out genome_filtered
plink --bfile genome_filtered --mind 0.1 --make-bed --out mind_filtered
plink --bfile mind_filtered --geno 0.02 --make-bed --out geno_filtered
plink --bfile geno_filtered --maf 0.05 --make-bed --out maf_filtered
plink --bfile maf_filtered --hwe 0.05 --make-bed --out qc_updated_phenotype_sex_all

# pca
plink --bfile merged_data --pca 10 --out pca_results
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' pca_results.eigenvec > pca_covar.txt

# logistic
plink --bfile qc_updated_phenotype_sex_all --logistic --allow-no-sex --out gwas_logistic_results_all
