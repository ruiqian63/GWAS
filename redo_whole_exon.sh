1. remove chr prefix in the chromosome of case cohort
sed 's/^#.*\n//; s/^chr//g' no_multi_case.vcf > no_multi_case_no_chr.vcf
sed 's/^##contig=<ID=chr/\n##contig=<ID=/g' no_multi_case_no_chr.vcf > no_multi_case_no_chr_cleaned.vcf
sed '/^$/d' no_multi_case_no_chr_cleaned.vcf > no_multi_case_no_chr_cleaned_no_empty.vcf
#result file: no_multi_case_no_chr_cleaned_no_empty.vcf.gz

2. filter snps (remove indels) in both cohort
bcftools view --types snps no_multi_case_no_chr_cleaned_no_empty.vcf.gz -Oz -o case_snps_only.vcf.gz
bcftools view --types snps pass_sfari_sibs_hg19toHg38.vcf.gz -Oz -o control_snps_only.vcf.gz

3. find shared snps in both (task1)
bcftools query -f '%CHROM\t%POS\t%ID\n' case_snps_only.vcf.gz | sort -u > case_snps.txt
#number: 7962572 
bcftools query -f '%CHROM\t%POS\t%ID\n' control_snps_only.vcf.gz | sort -u > control_snps.txt
#number: 3353360
comm -12 case_snps.txt control_snps.txt > shared_snps.txt
#number：2229473
#remove .
awk '$3 != "." {print $3}' shared_snps.txt > shared_snps_rsid.txt
#number：1925667

4. filter shared snp (task2)
bcftools index case_snps_only.vcf.gz
bcftools index control_snps_only.vcf.gz
bcftools view -R shared_snps_position.txt case_snps_only.vcf.gz -Oz -o case_snps_only_shared.vcf.gz
bcftools view -R shared_snps_position.txt control_snps_only.vcf.gz -Oz -o control_snps_only_shared.vcf.gz

5. change vcf header (task3)
bcftools view -h case_snps_only_shared.vcf.gz > hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt case_snps_only_shared.vcf.gz > case_snps_only_shared_removed_RG.vcf.gz
echo "removed_RG in case"

# change vcf header
bcftools view -h control_snps_only_shared.vcf.gz > hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt control_snps_only_shared.vcf.gz > control_snps_only_shared_removed_RG.vcf.gz
echo "removed_RG in control"

6. merge
bcftools index case_snps_only_shared_removed_RG.vcf.gz
#number of SNPs:	2071037
bcftools index control_snps_only_shared_removed_RG.vcf.gz
#number of SNPs:	1991100
bcftools merge case_snps_only_shared_removed_RG.vcf.gz control_snps_only_shared_removed_RG.vcf.gz -O z -o merged_1.vcf.gz
#number of SNPs:	2087684

7. add phenotype information
awk '{print $2, 2}' case.fam > case_ids.txt
awk '{print $2, 1}' control.fam > control_ids.txt
cat case_ids.txt control_ids.txt > all_ids.txt
awk 'NR==FNR{a[$1]=$2; next} {if ($2 in a) $6=a[$2]; print}' all_ids.txt merged.fam > updated_phenotype.fam
mv merged.bed updated_phenotype.bed
mv merged.bim updated_phenotype.bim

8. quality control
#plink --bfile updated_phenotype --genome --min 0.1 --allow-no-sex --make-bed --out genome_filtered
#2087684 variants and 3869 people pass filters and QC.
#3156 are cases and 713 are controls.
#Excluding 41546 variants on non-autosomes from IBD calculation.
plink --bfile updated_phenotype --mind 0.1 --make-bed --out mind_filtered
#Total genotyping rate in remaining samples is 0.978783.
#2087684 variants and 3868 people pass filters and QC.
#Among remaining phenotypes, 3156 are cases and 712 are controls.
plink --bfile mind_filtered --geno 0.1 --make-bed --out geno_filtered
#Total genotyping rate is 0.978783.
#16771 variants removed due to missing genotype data (--geno).
#2070913 variants and 3868 people pass filters and QC.
#Among remaining phenotypes, 3156 are cases and 712 are controls.
plink --bfile geno_filtered --maf 0.05 --make-bed --out maf_filtered
#Total genotyping rate is 0.98518.
#1984071 variants removed due to minor allele threshold(s)
#86842 variants and 3868 people pass filters and QC.
#Among remaining phenotypes, 3156 are cases and 712 are controls.
plink --bfile maf_filtered --hwe 0.05 --make-bed --out qc_updated_phenotype_all
#Total genotyping rate is 0.999267.
#16776 variants removed due to Hardy-Weinberg exact test.
#60029 variants and 3868 people pass filters and QC.
#Among remaining phenotypes, 3156 are cases and 712 are controls.

for EUR population and 1e-5 hwe:
75066 variants and 2012 people pass filters and QC.
Among remaining phenotypes, 1632 are cases and 380 are controls.

$remove chr6/x
plink --bfile qc_updated_phenotype_all --chr 1-5 7-22 --make-bed --out qc_updated_phenotype_all_no_chr6X


9. pca
plink --bfile qc_updated_phenotype_all --pca 10 --out pca_results
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' pca_results.eigenvec > pca_covar_10.txt

plink --bfile qc_updated_phenotype_all_no_chr6X --pca 10 --out pca_results_no_chr6X
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' pca_results_no_chr6X.eigenvec > pca_covar_10_no_chr6X.txt

plink --bfile qc_updated_phenotype_all --pca 3 --out pca_results
awk '{print $1, $2, $3, $4, $5}' pca_results.eigenvec > pca_covar_3.txt


10. logistic
plink --bfile qc_updated_phenotype_all --logistic --allow-no-sex --out gwas_logistic_results_all --adjust
plink --bfile qc_updated_phenotype_all --logistic --covar pca_covar_10.txt --allow-no-sex --out gwas_logistic_results_all_pca --adjust
plink --bfile qc_updated_phenotype_all --logistic --covar pca_covar_3.txt --allow-no-sex --out gwas_logistic_results_all_pca3 --adjust
plink2 --bfile qc_updated_phenotype_all --glm --covar pca_covar_10.txt --covar-variance-standardize --allow-no-sex --out gwas_logistic_results_all_pca10_plink2 

plink2 --bfile qc_updated_phenotype_all_no_chr6X --glm --covar pca_covar_10_no_chr6X.txt --covar-variance-standardize --allow-no-sex --out gwas_whole_exon_pca3_no_chr6X
awk 'NR==1 || $11 == "ADD" {print $0}' gwas_whole_exon_pca3_no_chr6X.PHENO1.glm.logistic.hybrid > filtered_gwas_whole_exon_pca3_no_chr6X.glm.logistic.hybrid


plink2 --bfile qc_updated_phenotype_all --glm --covar pca_covar_3.txt --covar-variance-standardize --allow-no-sex --out gwas_logistic_results_all_pca3_plink2 
awk 'NR==1 || $11 == "ADD" {print $0}' gwas_logistic_results_all_pca3_plink2.PHENO1.glm.logistic.hybrid > filtered_results.glm.logistic.hybrid


scp -r ch258782@e3-login.tch.harvard.edu:/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/gwas_whole_parameter/filtered_gwas_whole_exon_pca3_no_chr6X.glm.logistic.hybrid /Users/Rui/documents/GWAS/redo_no_chr6X/all 
