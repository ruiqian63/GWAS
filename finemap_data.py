import gwaslab as gl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sumstats = gl.Sumstats("/Users/Rui/documents/GWAS/redo/finemap_folate/filtered_results.glm.logistic.hybrid",fmt="plink2")
#print(sumstats.df)  # 打印数据表格

sumstats.basic_check()
sumstats.get_lead()
sumstats.plot_mqq()
plt.savefig('/Users/Rui/documents/GWAS/redo/finemap_folate_whole/mqq_plot.png', dpi=300)

locus = sumstats.filter_value('CHR==6 & POS>31000000 & POS<31500000')
locus.fill_data(to_fill=["BETA"])
locus.data
locus.harmonize(basic_check=False, ref_seq="/Users/Rui/documents/GWAS/redo/finemap_folate/human_g1k_v37.fasta")
locus.data
locus.data.to_csv("sig_locus.tsv",sep="\t",index=None)
locus.data["SNPID"].to_csv("sig_locus.snplist",sep="\t",index=None,header=None)
