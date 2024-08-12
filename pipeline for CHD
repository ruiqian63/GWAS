# merge vcf
bcftools index qc_control.vcf.gz
bcftools index qc_case.vcf.gz
bcftools merge qc_control.vcf.gz qc_case.vcf.gz -O z -o merged_filtered.vcf.gz

# change vcf header
bcftools view -h folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz > hdr.txt
# use sed to change header 
sed -i 's/ID=AD,Number=R,/ID=AD,Number=.,/' hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt folate_pass_coding_sfari_sibs_hg19toHg38.vcf.gz > folate_pass_coding_sfari_sibs_removed_RG.vcf.gz

# convert .vcf files
plink --vcf input_file.vcf.gz --make-bed --out output_prefix

# add sex information
awk '{if(NR>1) print $1, $2, $4}' file.sexcheck > updated_sex.txt
plink --bfile original_data --update-sex updated_sex.txt --make-bed --out updated_data

# add phenotype CHD
awk '{print $2, 2}' qc_case.fam > case_ids.txt
awk '{print $2, 1}' qc_control.fam > control_ids.txt
cat case_ids.txt control_ids.txt > all_ids.txt
awk 'NR==FNR{a[$1]=$2; next} {if ($2 in a) $6=a[$2]; print}' all_ids.txt updated_sex_all.fam > updated_phenotype.fam

# quality control
plink --bfile updated_phenotype_sex_all --genome --min 0.2 --allow-no-sex --make-bed --out genome_filtered
plink --bfile genome_filtered --mind 0.02 --make-bed --out mind_filtered
plink --bfile mind_filtered --geno 0.02 --make-bed --out geno_filtered
plink --bfile geno_filtered --maf 0.05 --make-bed --out maf_filtered
plink --bfile maf_filtered --hwe 1e-6 --make-bed --out qc_updated_phenotype_sex_all

# pca
plink --bfile merged_data --pca 10 --out pca_results
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' pca_results.eigenvec > pca_covar.txt

# logistic
plink --bfile qc_data --logistic --covar pca_covar.txt --pheno converted_phenotype.txt --allow-no-sex --out gwas_logistic_results
