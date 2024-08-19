#!/bin/bash

#split txt file
awk '$3 == "TRUE"' folate_cohort.txt > post_folate.txt
awk '$3 == "FALSE"' folate_cohort.txt > pre_folate.txt
cut -f1 post_folate.txt > post_folate_iids.txt
cut -f1 pre_folate.txt > pre_folate_iids.txt

#filter the merged cohort by folate
bcftools view -S post_folate_iids.txt -Oz -o post_folate.vcf.gz merged_all.vcf.gz
bcftools view -S pre_folate_iids.txt -Oz -o pre_folate.vcf.gz merged_all.vcf.gz

#construct phenotype_list
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' pre_folate.txt > phenotype_map.txt
awk '{if($2 == "FALSE") print $1, 1; else if($2 == "TRUE") print $1, 2}' post_folate.txt > phenotype_map.txt

#add phenotype
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' phenotype_map.txt pre_folate.fam > updated_phenotype.fam
awk 'NR==FNR {pheno[$1]=$2; next} {if ($1 in pheno) $6=pheno[$1]; print}' phenotype_map.txt post_folate.fam > updated_phenotype.fam
