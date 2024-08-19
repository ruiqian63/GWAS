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

#add sex
awk '{if(NR>1) print $1, $2, $4}' pcgc-exomes-n20631.vt.m1alt.U3.vqsrsi.fxh.anno.snpeff.dbnsfp.ESM1b.pf.plink.sexcheck > updated_sex1.txt
awk '{if(NR>1) print $1, $2, $4}' WES3-20160826.targets.a.f-016.m.b.f.snpeff.dbnsfp.1kg.esp.exac.trio.snpeff.plink.sexcheck > updated_sex.txt
plink --bfile updated_phenotype --update-sex updated_sex1.txt --make-bed --out updated_sex1
plink --bfile updated_sex1 --update-sex updated_sex.txt --make-bed --out updated_phenotype_sex

#quality contory
plink --bfile updated_phenotype_sex --genome --min 0.1 --allow-no-sex --make-bed --out genome_filtered
plink --bfile genome_filtered --mind 0.1 --make-bed --out mind_filtered
plink --bfile mind_filtered --geno 0.02 --make-bed --out geno_filtered
plink --bfile geno_filtered --maf 0.05 --make-bed --out maf_filtered
plink --bfile maf_filtered --hwe 0.05 --make-bed --out qc_updated_phenotype_sex

#pca

#GWAS
plink --bfile qc_updated_phenotype_sex --logistic --allow-no-sex --out gwas_logistic_results_pre_folate
