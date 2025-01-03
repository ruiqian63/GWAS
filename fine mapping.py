import gwaslab as gl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()

# 读取输入文件
df = pd.read_csv("/Users/Rui/documents/GWAS/redo/finemap/sig_locus.tsv", sep="\t")
ld = pd.read_csv("/Users/Rui/documents/GWAS/redo/finemap/sig_locus_mt.ld", sep="\t", header=None)
ld2 = pd.read_csv("/Users/Rui/documents/GWAS/redo/finemap/sig_locus_mt_r2.ld", sep="\t", header=None)

# 转换成 numpy 格式
R_df = ld.values
R_df2 = ld2.values

# 绘制 LD 矩阵的热图
plt.figure(figsize=(10,10), dpi=200)
fig, ax = plt.subplots(ncols=2, figsize=(20,10))
sns.heatmap(data=R_df, cmap="Spectral", ax=ax[0])
sns.heatmap(data=R_df2, ax=ax[1])
ax[0].set_title("LD r matrix")
ax[1].set_title("LD r2 matrix")

# 保存热图
plt.savefig("/Users/Rui/documents/GWAS/redo/finemap/heatmap_ld.png", dpi=200)

# 继续 SuSiE 分析
susieR = importr('susieR')
ro.r('set.seed(123)')
fit = susieR.susie_rss(
    bhat = df["BETA"].values.reshape((len(R_df), 1)),
    shat = df["SE"].values.reshape((len(R_df), 1)),
    R = R_df,
    L = 10,
    n = 503
)

# 获取 CS 信息并更新 DataFrame
df["cs"] = 0
n_cs = len(susieR.susie_get_cs(fit, coverage=0.95, min_abs_corr=0.5, Xcorr=R_df)[0])
for i in range(n_cs):
    cs_index = susieR.susie_get_cs(fit, coverage=0.95, min_abs_corr=0.5, Xcorr=R_df)[0][i]
    df.loc[np.array(cs_index)-1, "cs"] = i + 1
df["pip"] = np.array(susieR.susie_get_pip(fit))

# 绘制 P 和 PIP 图
fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(15, 7), height_ratios=(4, 1))
df["MLOG10P"] = -np.log10(df["P"])
col_to_plot = "MLOG10P"
p = axes[0].scatter(df["POS"], df[col_to_plot], c=ld[df["P"].idxmin()]**2)

axes[0].scatter(df.loc[df["cs"] == 1, "POS"], df.loc[df["cs"] == 1, col_to_plot],
                marker='o', s=40, c="None", edgecolors='black', label="Variants in credible set 1")
axes[0].scatter(df.loc[(df["CHR"] == 2) & (df["POS"] == 55620927), "POS"],
                df.loc[(df["CHR"] == 2) & (df["POS"] == 55620927), col_to_plot],
                marker='x', s=40, c="red", edgecolors='black', label="Causal")

plt.colorbar(p, label="Rsq with the lead variant")
axes[0].set_xlabel("position")
axes[0].set_xlim((30269400, 33269400))
axes[0].set_ylabel(col_to_plot)
axes[0].legend()

p = axes[1].scatter(df["POS"], df["pip"], c=ld[df["P"].idxmin()]**2)
axes[1].scatter(df.loc[df["cs"] == 1, "POS"], df.loc[df["cs"] == 1, "pip"],
                marker='o', s=40, c="None", edgecolors='black', label="Variants in credible set 1")

plt.colorbar(p, label="Rsq with the lead variant")
axes[1].set_xlabel("position")
axes[1].set_xlim((30269400, 33269400))
axes[1].set_ylabel("PIP")
axes[1].legend()

# 保存 P 和 PIP 图像
plt.savefig("/Users/Rui/documents/GWAS/redo/finemap/p_pip_plot.png", dpi=200)

plt.show()
