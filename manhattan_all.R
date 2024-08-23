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

# 确保染色体是因子，以便保持顺序
data$CHR <- factor(data$CHR, levels = unique(data$CHR))

# 为每个染色体分配一个基因组位置范围
chromosome_info <- data %>%
  group_by(CHR) %>%
  summarize(chr_len = max(BP)) %>%
  mutate(tot = cumsum(chr_len) - chr_len)

# 计算缩放因子（缩小前几个染色体的范围）
scale_factor <- 0.5 # 例如，将前几个染色体的范围缩小为原来的50%
n_small_chromosomes <- 5 # 假设前5个染色体需要缩小

data <- data %>%
  inner_join(chromosome_info, by = "CHR") %>%
  mutate(BP_cum = ifelse(as.numeric(CHR) <= n_small_chromosomes,
                         (BP / max(BP)) * scale_factor * max(chr_len) + tot,
                         BP + tot))

# 创建染色体标签的位置
axis_set <- chromosome_info %>%
  mutate(center = tot + chr_len / 2,
         scaled_center = ifelse(as.numeric(CHR) <= n_small_chromosomes,
                                (chr_len / max(chr_len)) * scale_factor * max(chr_len) / 2 + tot,
                                center))

# 设置颜色
colors <- rep(c("blue", "orange"), length.out = length(unique(data$CHR)))

# 创建图形
ggplot(data, aes(x = BP_cum, y = `-log10(p)`, color = as.factor(CHR))) +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  labs(x = "Genomic Position", y = "-log10(p-value)", title = "Manhattan Plot") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = axis_set$scaled_center, labels = axis_set$CHR)
