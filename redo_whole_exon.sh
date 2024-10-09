#1. remove chr prefix in the chromosome of case cohort
sed 's/^#.*\n//; s/^chr//g' no_multi_case.vcf > no_multi_case_no_chr.vcf
sed 's/^##contig=<ID=chr/\n##contig=<ID=/g' no_multi_case_no_chr.vcf > no_multi_case_no_chr_cleaned.vcf
sed '/^$/d' no_multi_case_no_chr_cleaned.vcf > no_multi_case_no_chr_cleaned_no_empty.vcf
#result file: no_multi_case_no_chr_cleaned_no_empty.vcf.gz

#2. filter snps (remove indels) in both cohort
bcftools view --types snps no_multi_case_no_chr_cleaned_no_empty.vcf.gz -Oz -o case_snps_only.vcf.gz
bcftools view --types snps pass_sfari_sibs_hg19toHg38.vcf.gz -Oz -o control_snps_only.vcf.gz

#3. find shared snps in both
bcftools query -f '%CHROM\t%POS\t%ID\n' case_snps_only.vcf.gz | sort -u > case_snps.txt
#number: 
bcftools query -f '%CHROM\t%POS\t%ID\n' control_snps_only.vcf.gz | sort -u > control_snps.txt
#number: 
comm -12 case_snps.txt control_snps.txt > shared_snps.txt
#number：

