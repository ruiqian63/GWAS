#1. remove chr prefix in the chromosome of case cohort
sed 's/^#.*\n//; s/^chr//g' no_multi_case.vcf > no_multi_case_no_chr.vcf
sed 's/^##contig=<ID=chr/\n##contig=<ID=/g' no_multi_case_no_chr.vcf > no_multi_case_no_chr_cleaned.vcf
sed '/^$/d' no_multi_case_no_chr_cleaned.vcf > no_multi_case_no_chr_cleaned_no_empty.vcf
#result file: no_multi_case_no_chr_cleaned_no_empty.vcf.gz

#2. filter snps (remove indels) in both cohort
