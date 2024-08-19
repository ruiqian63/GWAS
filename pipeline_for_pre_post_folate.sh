#!/bin/bash

#split txt file
awk '$3 == "TRUE"' folate_cohort.txt > post_folate.txt
awk '$3 == "FALSE"' folate_cohort.txt > pre_folate.txt
cut -f1 post_folate.txt > post_folate_iids.txt
cut -f1 pre_folate.txt > pre_folate_iids.txt

#filter the merged cohort by folate
bcftools view -S post_folate_iids.txt -Oz -o post_folate.vcf.gz merged_all.vcf.gz
bcftools view -S pre_folate_iids.txt -Oz -o pre_folate.vcf.gz merged_all.vcf.gz
