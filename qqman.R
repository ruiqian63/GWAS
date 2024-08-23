# 加载qqman包
library(qqman)

# 读取GWAS结果
gwasResults <- read.table("E:\\BCH_research\\GWAS_exp\\gwas_logistic_results.assoc.logistic", header=TRUE)

# 绘制曼哈顿图
manhattan(gwasResults, main="Manhattan Plot", ylim=c(0, -log10(1e-8)))

# 绘制QQ图
qq(gwasResults$P, main="Q-Q plot of GWAS p-values")
