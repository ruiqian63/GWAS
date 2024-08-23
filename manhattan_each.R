library(ggplot2)
library(dplyr)

# 读取数据文件
file_path <- "E:/BCH_research/GWAS_exp/gwas_logistic_results.assoc.logistic"
data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)

# 假设数据文件包含以下列：CHR, BP, P
# CHR: 染色体编号
# BP: 基因组位置
# P: p值

# 转换p值为-log10(p)
data <- data %>%
  mutate(`-log10(p)` = -log10(P))

# 设置颜色
colors <- rep(c("blue", "orange"), length.out = length(unique(data$CHR)))

# 创建图形
ggplot(data, aes(x = BP, y = `-log10(p)`, color = as.factor(CHR))) +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  labs(x = "Genomic Position", y = "-log10(p-value)", title = "Manhattan Plot") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~CHR, scales = "free_x")
