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
bcftools index control_snps_only_shared_removed_RG.vcf.gz
bcftools merge case_snps_only_shared_removed_RG.vcf.gz control_snps_only_shared_removed_RG.vcf.gz -O z -o merged_1.vcf.gz
