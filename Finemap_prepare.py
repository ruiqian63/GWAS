import gwaslab as gl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

awk 'NR==1 || $11 == "ADD" {print $0}' post_gwas_logistic_results_all_pca10_plink2.PHENO1.glm.logistic.hybrid > post_filtered_results.glm.logistic.hybrid

sumstats = gl.Sumstats("./pre_gwas_logistic_results_all_pca10_plink2.PHENO1.glm.logistic.hybrid",fmt="plink2")
sumstats.basic_check()
sumstats.get_lead()
sumstats.plot_mqq()
plt.savefig('mqq_plot.png')

locus = sumstats.filter_value('CHR==6 & POS>30269400 & POS<33269400')
locus.fill_data(to_fill=["BETA"])
locus.data
locus.harmonize(basic_check=False, ref_seq="./human_g1k_v37.fasta")
locus.data
locus.data.to_csv("sig_locus.tsv",sep="\t",index=None)
locus.data["SNPID"].to_csv("sig_locus.snplist",sep="\t",index=None,header=None)

plink \
  --bfile "./post_qc_updated_phenotype_all" \
  --keep-allele-order \
  --r square \
  --extract sig_locus.snplist \
  --out sig_locus_mt

plink \
  --bfile "./post_qc_updated_phenotype_all" \
  --keep-allele-order \
  --r2 square \
  --extract sig_locus.snplist \
  --out sig_locus_mt_r2
