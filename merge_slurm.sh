#!/bin/bash

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=48:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=liftover # Job name
#SBATCH --output=hostname.%j.out # Name of the output file
#SBATCH --error=hostname.%j.error
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of threads/tasks on one node
#SBATCH --mem=12G

## commands
source /programs/biogrids.shrc

#get vcf index
bcftools index no_multi_case.vcf.gz
bcftools index no_multi_control.vcf.gz

# finding shared variants
bcftools query -f '%CHROM\t%POS\n' no_multi_control.vcf.gz | sort -u > file1_positions.txt
bcftools query -f '%CHROM\t%POS\n' no_multi_case.vcf.gz | sort -u > file2_positions.txt
comm -12 <(sort file1_positions.txt) <(sort file2_positions.txt) > shared_positions.txt

# filter VCF files based on shared positions
vcftools --gzvcf no_multi_control.vcf.gz --positions shared_positions.txt --recode --stdout | bgzip -c > filtered_no_multi_control.vcf.gz
vcftools --gzvcf no_multi_case.vcf.gz --positions shared_positions.txt --recode --stdout | bgzip -c > filtered_no_multi_case.vcf.gz

#get vcf index
bcftools index filtered_no_multi_case.vcf.gz
bcftools index filtered_no_multi_control.vcf.gz

# change vcf header
bcftools view -h filtered_no_multi_case.vcf.gz > hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt filtered_no_multi_case.vcf.gz > filtered_no_multi_case_removed_RG.vcf.gz

# change vcf header
bcftools view -h filtered_no_multi_control.vcf.gz > hdr.txt
# change to all instances of Number=R/G in the header
sed -i 's/Number=R/Number=./g' hdr.txt
sed -i 's/Number=G/Number=./g' hdr.txt
# use the modified header to update the original VCF file
bcftools reheader -h hdr.txt filtered_no_multi_control.vcf.gz > filtered_no_multi_control_removed_RG.vcf.gz

#get vcf index
bcftools index filtered_no_multi_case_removed_RG.vcf.gz
bcftools index filtered_no_multi_control_removed_RG.vcf.gz

# merge vcf
bcftools merge filtered_no_multi_case_removed_RG.vcf.gz filtered_no_multi_control_removed_RG.vcf.gz -O z -o merged_whole.vcf.gz
