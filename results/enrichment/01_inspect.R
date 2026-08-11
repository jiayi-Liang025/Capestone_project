# ============================================================
# Step 1
# 检查 6 个 CSV 文件的结构
# ============================================================

# ------------------------------------------------------------
# 1. 设置工作目录
# ------------------------------------------------------------

# 方法 A：如果你已经把 RStudio Project 建好了
# 可以直接使用：
getwd()

# 如果 getwd() 不是你的项目目录，
# 可以用下面这种方式修改：
#
# setwd("你的/CMap_drug_repurposing/路径")


# ------------------------------------------------------------
# 2. 找到 6 个 CSV
# ------------------------------------------------------------

raw_dir <- "results/enrichment"

files <- list.files(
  path = raw_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

files

# ============================================================
# Step 2
# 从 GO enrichment 结果中提取 geneID
# ============================================================

raw_dir <- "results/enrichment"

files <- list.files(
  path = raw_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# 查看每个文件是否存在 geneID 列
for (f in files) {
  
  df <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("\n============================================\n")
  cat(basename(f), "\n")
  cat("============================================\n")
  
  cat("列名：\n")
  print(colnames(df))
  
  # 自动寻找 geneID
  gene_col <- colnames(df)[
    tolower(colnames(df)) %in% c(
      "geneid",
      "gene_id",
      "geneids",
      "gene"
    )
  ]
  
  if (length(gene_col) > 0) {
    cat("\n找到基因列：", gene_col[1], "\n")
    print(head(df[[gene_col[1]]], 5))
  } else {
    cat("\n⚠️ 没有找到 geneID 列\n")
  }
}

# ============================================================
# Step 2
# 从 GO enrichment 结果中提取唯一 Gene Symbol
# ============================================================

raw_dir <- "results/enrichment"

# 输出目录
gene_dir <- "results/gene_sets"

if (!dir.exists(gene_dir)) {
  dir.create(gene_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. 找到 6 个 GO enrichment 文件
# ------------------------------------------------------------

files <- list.files(
  path = raw_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

print(basename(files))


# ------------------------------------------------------------
# 2. 提取 geneID
# ------------------------------------------------------------

extract_genes <- function(file) {
  
  df <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # 检查 geneID
  if (!"geneID" %in% colnames(df)) {
    stop(
      paste("没有找到 geneID:", basename(file))
    )
  }
  
  # 取 geneID
  gene_strings <- df$geneID
  
  # 去除 NA
  gene_strings <- gene_strings[
    !is.na(gene_strings)
  ]
  
  # 合并成一个字符串
  gene_strings <- paste(
    gene_strings,
    collapse = "/"
  )
  
  # 按 / 拆分
  genes <- unlist(
    strsplit(
      gene_strings,
      split = "/"
    )
  )
  
  # 去除空格
  genes <- trimws(genes)
  
  # 去除空值
  genes <- genes[
    genes != "" &
      !is.na(genes)
  ]
  
  # 去重
  genes <- unique(genes)
  
  # 排序
  genes <- sort(genes)
  
  return(genes)
}


# ------------------------------------------------------------
# 3. 对 6 个文件分别提取
# ------------------------------------------------------------

gene_sets <- list()

for (file in files) {
  
  file_name <- tools::file_path_sans_ext(
    basename(file)
  )
  
  cat("\n----------------------------------------\n")
  cat(file_name, "\n")
  
  genes <- extract_genes(file)
  
  gene_sets[[file_name]] <- genes
  
  cat("Unique genes:", length(genes), "\n")
}


# ------------------------------------------------------------
# 4. 保存结果
# ------------------------------------------------------------

for (name in names(gene_sets)) {
  
  output_file <- file.path(
    gene_dir,
    paste0(name, "_genes.csv")
  )
  
  write.csv(
    data.frame(
      Gene = gene_sets[[name]]
    ),
    output_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  cat(
    "Saved:",
    output_file,
    "\n"
  )
}


# ------------------------------------------------------------
# 5. 汇总
# ------------------------------------------------------------

gene_summary <- data.frame(
  File = names(gene_sets),
  Unique_Genes = sapply(
    gene_sets,
    length
  )
)

print(gene_summary)


# 保存汇总
write.csv(
  gene_summary,
  "results/gene_sets/gene_set_summary.csv",
  row.names = FALSE
)

# 查看整个 results 目录下面有哪些 CSV

all_csv <- list.files(
  "results",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

print(all_csv)
gene_summary_check <- read.csv(
  "results/gene_sets/gene_set_summary.csv",
  stringsAsFactors = FALSE
)

cat("行数:", nrow(gene_summary_check), "\n")
cat("列名:\n")
print(names(gene_summary_check))

cat("\n前几行:\n")
print(head(gene_summary_check))

gene_summary_check

#
list.files(
  "results/gene_sets",
  full.names = TRUE
)
#
gene_test <- read.csv(
  "results/gene_sets/COVID19_vs_Healthy_GO_down_genes.csv",
  stringsAsFactors = FALSE
)

print(names(gene_test))
print(head(gene_test))
#
# 6 个 gene set 文件
gene_files <- c(
  "COVID19_vs_Healthy_GO_down",
  "COVID19_vs_Healthy_GO_up",
  "COVID19_vs_OtherVirus_GO_down",
  "COVID19_vs_OtherVirus_GO_up",
  "OtherVirus_vs_Healthy_GO_down",
  "OtherVirus_vs_Healthy_GO_up"
)

# 读取成 gene sets
gene_sets <- lapply(gene_files, function(x) {
  dat <- read.csv(
    paste0("results/gene_sets/", x, "_genes.csv"),
    stringsAsFactors = FALSE
  )
  
  unique(dat$Gene)
})

# 给每个 gene set 命名
names(gene_sets) <- gene_files

# 查看每组基因数量
sapply(gene_sets, length)
#

overlap_matrix <- outer(
  names(gene_sets),
  names(gene_sets),
  Vectorize(function(a, b) {
    length(intersect(gene_sets[[a]], gene_sets[[b]]))
  })
)

rownames(overlap_matrix) <- names(gene_sets)
colnames(overlap_matrix) <- names(gene_sets)

print(overlap_matrix)
#
common_genes <- Reduce(intersect, gene_sets)

cat("6组共同基因数量:", length(common_genes), "\n")

print(common_genes)
#

genes_95 <- intersect(
  gene_sets[["COVID19_vs_Healthy_GO_up"]],
  gene_sets[["COVID19_vs_OtherVirus_GO_down"]]
)

cat("共同基因数量:", length(genes_95), "\n")
print(genes_95)
#
# 1. 95 genes
genes_95 <- intersect(
  gene_sets[["COVID19_vs_Healthy_GO_up"]],
  gene_sets[["COVID19_vs_OtherVirus_GO_down"]]
)

# 2. 60 genes
genes_60 <- intersect(
  gene_sets[["COVID19_vs_Healthy_GO_up"]],
  gene_sets[["OtherVirus_vs_Healthy_GO_up"]]
)

# 3. 51 genes
genes_51 <- intersect(
  gene_sets[["COVID19_vs_OtherVirus_GO_down"]],
  gene_sets[["OtherVirus_vs_Healthy_GO_up"]]
)

# 4. 6 genes
genes_6 <- intersect(
  gene_sets[["COVID19_vs_OtherVirus_GO_up"]],
  gene_sets[["OtherVirus_vs_Healthy_GO_down"]]
)

cat("95 genes:", length(genes_95), "\n")
cat("60 genes:", length(genes_60), "\n")
cat("51 genes:", length(genes_51), "\n")
cat("6 genes:", length(genes_6), "\n")
#
cat("===== 60 genes =====\n")
print(genes_60)

cat("\n===== 51 genes =====\n")
print(genes_51)

cat("\n===== 6 genes =====\n")
print(genes_6)

#visual
# ================================
# Jaccard similarity
# ================================

jaccard_matrix <- matrix(
  0,
  nrow = length(gene_sets),
  ncol = length(gene_sets)
)

rownames(jaccard_matrix) <- names(gene_sets)
colnames(jaccard_matrix) <- names(gene_sets)

for (i in seq_along(gene_sets)) {
  for (j in seq_along(gene_sets)) {
    
    intersection <- length(
      intersect(gene_sets[[i]], gene_sets[[j]])
    )
    
    union_genes <- length(
      union(gene_sets[[i]], gene_sets[[j]])
    )
    
    jaccard_matrix[i, j] <- intersection / union_genes
  }
}

# 查看 Jaccard matrix
print(round(jaccard_matrix, 3))


#
# 在整个项目中搜索包含 all_DEGs 的文件
all_csv_files <- list.files(
  path = ".",
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$"
)

print(all_csv_files)

#
# ============================================
# Read the three DEG results
# ============================================

deg_healthy <- read.csv(
  "results/differential_expression/COVID19_vs_Healthy_all_genes.csv",
  stringsAsFactors = FALSE
)

deg_othervirus <- read.csv(
  "results/differential_expression/COVID19_vs_OtherVirus_all_genes.csv",
  stringsAsFactors = FALSE
)

deg_ov_healthy <- read.csv(
  "results/differential_expression/OtherVirus_vs_Healthy_all_genes.csv",
  stringsAsFactors = FALSE
)

# Check column names
names(deg_healthy)
names(deg_othervirus)
names(deg_ov_healthy)
#

# ============================================
# Prepare log2FC matrix for correlation
# ============================================

# 只保留 gene 和 log2FoldChange
fc_healthy <- deg_healthy[, c("gene", "log2FoldChange")]
fc_othervirus <- deg_othervirus[, c("gene", "log2FoldChange")]
fc_ov_healthy <- deg_ov_healthy[, c("gene", "log2FoldChange")]

# 改列名
colnames(fc_healthy)[2] <- "COVID19_vs_Healthy"
colnames(fc_othervirus)[2] <- "COVID19_vs_OtherVirus"
colnames(fc_ov_healthy)[2] <- "OtherVirus_vs_Healthy"

# 合并三个 comparison
fc_all <- merge(
  fc_healthy,
  fc_othervirus,
  by = "gene"
)

fc_all <- merge(
  fc_all,
  fc_ov_healthy,
  by = "gene"
)

# 看结果
cat("共同基因数量:", nrow(fc_all), "\n")
head(fc_all)
# ============================================
# Correlation of gene expression signatures
# ============================================

fc_matrix <- fc_all[, c(
  "COVID19_vs_Healthy",
  "COVID19_vs_OtherVirus",
  "OtherVirus_vs_Healthy"
)]

# Pearson correlation
cor_pearson <- cor(
  fc_matrix,
  method = "pearson",
  use = "pairwise.complete.obs"
)

cat("===== Pearson correlation =====\n")
print(round(cor_pearson, 3))


# Spearman correlation
cor_spearman <- cor(
  fc_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

cat("\n===== Spearman correlation =====\n")
print(round(cor_spearman, 3))
cat("===== Spearman correlation =====\n")
print(round(cor_spearman, 3))


#

# ============================================
# Correlation heatmap
# ============================================

library(ggplot2)

cor_df <- as.data.frame(as.table(cor_spearman))

colnames(cor_df) <- c(
  "Comparison1",
  "Comparison2",
  "Correlation"
)

ggplot(
  cor_df,
  aes(
    x = Comparison1,
    y = Comparison2,
    fill = Correlation
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", Correlation)),
    size = 5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Spearman Correlation of Gene Expression Signatures",
    x = NULL,
    y = NULL,
    fill = "Spearman r"
  )





# ============================================
# STEP 1: DEG overlap
# ============================================

# 取显著 DEGs
deg_COVID_H <- deg_healthy$gene[!is.na(deg_healthy$padj) & deg_healthy$padj < 0.05]

deg_COVID_OV <- deg_othervirus$gene[
  !is.na(deg_othervirus$padj) & deg_othervirus$padj < 0.05
]

deg_OV_H <- deg_ov_healthy$gene[
  !is.na(deg_ov_healthy$padj) & deg_ov_healthy$padj < 0.05
]

# 去除重复
deg_COVID_H <- unique(deg_COVID_H)
deg_COVID_OV <- unique(deg_COVID_OV)
deg_OV_H <- unique(deg_OV_H)

# 看每组有多少
cat("COVID19 vs Healthy:", length(deg_COVID_H), "\n")
cat("COVID19 vs OtherVirus:", length(deg_COVID_OV), "\n")
cat("OtherVirus vs Healthy:", length(deg_OV_H), "\n")



# ============================================
# DEG overlap
# ============================================

overlap_COVID_H_COVID_OV <- intersect(
  deg_COVID_H,
  deg_COVID_OV
)

overlap_COVID_H_OV_H <- intersect(
  deg_COVID_H,
  deg_OV_H
)

overlap_COVID_OV_OV_H <- intersect(
  deg_COVID_OV,
  deg_OV_H
)

overlap_all_three <- Reduce(
  intersect,
  list(
    deg_COVID_H,
    deg_COVID_OV,
    deg_OV_H
  )
)

cat("\n===== DEG OVERLAP =====\n")

cat(
  "COVID19_vs_Healthy ∩ COVID19_vs_OtherVirus:",
  length(overlap_COVID_H_COVID_OV),
  "\n"
)

cat(
  "COVID19_vs_Healthy ∩ OtherVirus_vs_Healthy:",
  length(overlap_COVID_H_OV_H),
  "\n"
)

cat(
  "COVID19_vs_OtherVirus ∩ OtherVirus_vs_Healthy:",
  length(overlap_COVID_OV_OV_H),
  "\n"
)

cat(
  "ALL THREE:",
  length(overlap_all_three),
  "\n"
)
cat("COVID19 vs Healthy:", nrow(deg_healthy), "\n")
cat("COVID19 vs OtherVirus:", nrow(deg_othervirus), "\n")
cat("OtherVirus vs Healthy:", nrow(deg_ov_healthy), "\n")

covid_genes <- deg_healthy$gene
othervirus_healthy_genes <- deg_ov_healthy$gene

covid_specific_genes <- setdiff(
  covid_genes,
  othervirus_healthy_genes
)

cat("COVID-specific genes:", length(covid_specific_genes), "\n")


# ==============================
# STEP 1: Define DEG gene sets
# ==============================

genes_COVID_H <- unique(deg_healthy$gene)
genes_COVID_OV <- unique(deg_othervirus$gene)
genes_OV_H <- unique(deg_ov_healthy$gene)

cat("COVID19 vs Healthy:", length(genes_COVID_H), "\n")
cat("COVID19 vs OtherVirus:", length(genes_COVID_OV), "\n")
cat("OtherVirus vs Healthy:", length(genes_OV_H), "\n")

# Three-way shared genes
shared_all_three <- Reduce(
  intersect,
  list(genes_COVID_H, genes_COVID_OV, genes_OV_H)
)

cat("Three-way shared genes:", length(shared_all_three), "\n")

# ============================================
# CORRECT DEG SETS FOR OVERLAP
# ============================================

sig_COVID_H <- read.csv(
  "./results/differential_expression/COVID19_vs_Healthy_significant_genes.csv"
)

sig_COVID_OV <- read.csv(
  "./results/differential_expression/COVID19_vs_OtherVirus_significant_genes.csv"
)

sig_OV_H <- read.csv(
  "./results/differential_expression/OtherVirus_vs_Healthy_significant_genes.csv"
)

# Check columns
names(sig_COVID_H)
names(sig_COVID_OV)
names(sig_OV_H)

#
genes_COVID_H <- unique(sig_COVID_H$gene)
genes_COVID_OV <- unique(sig_COVID_OV$gene)
genes_OV_H <- unique(sig_OV_H$gene)

cat("Significant COVID19 vs Healthy:",
    length(genes_COVID_H), "\n")

cat("Significant COVID19 vs OtherVirus:",
    length(genes_COVID_OV), "\n")

cat("Significant OtherVirus vs Healthy:",
    length(genes_OV_H), "\n")
# Three-way overlap
shared_all_three <- Reduce(
  intersect,
  list(
    genes_COVID_H,
    genes_COVID_OV,
    genes_OV_H
  )
)

cat(
  "Three-way shared significant DEGs:",
  length(shared_all_three),
  "\n"
)

# ============================================
# DEG OVERLAP — UPSET PLOT
# ============================================

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}

library(UpSetR)

# Create named list of significant DEG sets

deg_sets <- list(
  COVID19_vs_Healthy = genes_COVID_H,
  COVID19_vs_OtherVirus = genes_COVID_OV,
  OtherVirus_vs_Healthy = genes_OV_H
)

# Convert to UpSet input
upset_data <- fromList(deg_sets)

# Check
head(upset_data)


png(
  "results/DEG_overlap_UpSet.png",
  width = 2400,
  height = 1800,
  res = 300
)

upset(
  upset_data,
  sets = c(
    "COVID19_vs_Healthy",
    "COVID19_vs_OtherVirus",
    "OtherVirus_vs_Healthy"
  ),
  order.by = "freq",
  mainbar.y.label = "Number of shared significant DEGs",
  sets.x.label = "Number of significant DEGs",
  text.scale = 1.5
)

dev.off()






length(shared_all_three)
head(shared_all_three)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)
#
gene_conversion <- bitr(
  shared_all_three,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_conversion)
#
cat("Original genes:", length(shared_all_three), "\n")
cat("Successfully mapped:", nrow(gene_conversion), "\n")

#
ego_three <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(as.data.frame(ego_three))
cat("Significant GO BP terms:", nrow(as.data.frame(ego_three)), "\n")
dotplot(
  ego_three,
  showCategory = 15,
  title = "GO Biological Process: Three-way Shared DEGs"
)
##Conclusion
##The genes consistently dysregulated across COVID-19, other viral infection, 
##and healthy comparisons are strongly enriched in antiviral innate immunity, 
##inflammatory signaling, interferon-related responses, and leukocyte recruitment.



covid_specific <- setdiff(
  genes_COVID_H,
  union(genes_COVID_OV, genes_OV_H)
)

cat("COVID19-specific significant DEGs:", length(covid_specific), "\n")
head(covid_specific)
cat("COVID19-specific significant DEGs:", length(covid_specific), "\n")


#
covid_conversion <- bitr(
  covid_specific,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("Original COVID19-specific genes:", length(covid_specific), "\n")
cat("Successfully mapped:", nrow(covid_conversion), "\n")

#GO enrichment
ego_covid_specific <- enrichGO(
  gene = covid_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
#numbers of pathways:Significant GO BP terms: 161 
cat(
  "Significant GO BP terms:",
  nrow(as.data.frame(ego_covid_specific)),
  "\n"
)

#
dotplot(
  ego_covid_specific,
  showCategory = 15,
  title = "GO Biological Process: COVID19-specific DEGs"
)


#
# ================================
# OtherVirus-specific significant DEGs
# ================================

other_specific <- setdiff(
  genes_OV_H,
  union(genes_COVID_H, genes_COVID_OV)
)

cat(
  "OtherVirus-specific significant DEGs:",
  length(other_specific),
  "\n"
)


# ==========================================
# GO enrichment: OtherVirus-specific DEGs
# ==========================================

ego_other <- enrichGO(
  gene          = other_specific,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 查看结果
head(ego_other)

# 保存 GO enrichment results
write.csv(
  as.data.frame(ego_other),
  "OtherVirus_specific_GO_BP.csv",
  row.names = FALSE
)

# ==========================================
# Plot
# ==========================================

library(stringr)

p_other_GO <- dotplot(
  ego_other,
  showCategory = 10
) +
  ggtitle("GO Biological Process: OtherVirus-specific DEGs") +
  scale_y_discrete(
    labels = function(x) {
      str_wrap(x, width = 32)
    }
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.y = element_text(
      size = 9,
      lineheight = 0.9
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.title = element_text(
      size = 11
    ),
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

print(p_other_GO)

ggsave(
  "OtherVirus_specific_GO_BP.png",
  p_other_GO,
  width = 12,
  height = 12,
  dpi = 300
)


##COVID19-specific vs OtherVirus-specific FUNCTION COMPARISON
library(dplyr)
library(ggplot2)
library(stringr)

# ================================
# COVID19-specific vs OtherVirus-specific GO
# ================================

# COVID19-specific：取 Top 10
covid_GO <- ego_covid_specific@result %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Group = "COVID19-specific"
  )

# OtherVirus-specific：取 Top 10
other_GO <- ego_other@result %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Group = "OtherVirus-specific"
  )

# 合并
GO_compare <- bind_rows(
  covid_GO,
  other_GO
)

# GeneRatio 转成 numeric
GO_compare <- GO_compare %>%
  mutate(
    GeneRatio_num = sapply(
      strsplit(as.character(GeneRatio), "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
  )

# 长 GO 名称自动换行
GO_compare <- GO_compare %>%
  mutate(
    Description = str_wrap(
      Description,
      width = 38
    )
  )

# 设置每组内部排序
GO_compare <- GO_compare %>%
  group_by(Group) %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(Description)
    )
  ) %>%
  ungroup()

# ================================
# Plot
# ================================

p_GO_compare <- ggplot(
  GO_compare,
  aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point() +
  facet_wrap(
    ~ Group,
    scales = "free_y"
  ) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    trans = "reverse"
  ) +
  labs(
    title = "GO Biological Process Comparison",
    x = "Gene Ratio",
    y = NULL,
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16
    ),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9,
      lineheight = 0.9
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.title.x = element_text(
      size = 12
    ),
    panel.grid.major.y = element_line(
      linewidth = 0.3
    ),
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

print(p_GO_compare)

# ================================
# 保存
# ================================

ggsave(
  "COVID19_vs_OtherVirus_GO_comparison.png",
  p_GO_compare,
  width = 15,
  height = 10,
  dpi = 300
)

#Test GO plot
library(enrichplot)
library(ggplot2)

p_covid <- dotplot(
  ego_covid_specific,
  showCategory = 10
) +
  ggtitle("COVID19-specific DEGs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

p_other <- dotplot(
  ego_other,
  showCategory = 10
) +
  ggtitle("OtherVirus-specific DEGs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

print(p_covid)
print(p_other)

#

library(dplyr)
library(ggplot2)
library(scales)

# ================================
# 1. 准备 COVID19-specific 数据
# ================================

covid_plot <- ego_covid_specific@result %>%
  slice_head(n = 10) %>%
  mutate(
    Group = "COVID19-specific",
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  ) %>%
  select(
    Group,
    Description,
    GeneRatio_num,
    Count,
    p.adjust
  )


# ================================
# 2. 准备 OtherVirus-specific 数据
# ================================

other_plot <- ego_other@result %>%
  slice_head(n = 10) %>%
  mutate(
    Group = "OtherVirus-specific",
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  ) %>%
  select(
    Group,
    Description,
    GeneRatio_num,
    Count,
    p.adjust
  )


# ================================
# 3. 合并
# ================================

go_compare <- bind_rows(
  covid_plot,
  other_plot
)


# 看一下数据是否正常
print(go_compare)


ego_covid_df <- ego_covid_specific@result
ego_other_df <- ego_other@result

ego_covid_df$GeneRatio_num <- sapply(
  strsplit(ego_covid_df$GeneRatio, "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)

ego_other_df$GeneRatio_num <- sapply(
  strsplit(ego_other_df$GeneRatio, "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)


#
library(patchwork)
covid_plot_df <- ego_covid_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = str_wrap(Description, width = 30),
    Description = factor(Description, levels = rev(Description))
  )

other_plot_df <- ego_other_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = str_wrap(Description, width = 30),
    Description = factor(Description, levels = rev(Description))
  )

#
p_covid <- ggplot(
  covid_plot_df,
  aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    trans = "reverse"
  ) +
  scale_size_continuous(
    range = c(3, 10)
  ) +
  labs(
    title = "COVID19-specific DEGs",
    x = "Gene Ratio",
    y = NULL,
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.text.y = element_text(
      size = 11
    ),
    axis.text.x = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    )
  )
print(p_covid)


#
p_other <- ggplot(
  other_plot_df,
  aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    trans = "reverse"
  ) +
  scale_size_continuous(
    range = c(3, 10)
  ) +
  labs(
    title = "OtherVirus-specific DEGs",
    x = "Gene Ratio",
    y = NULL,
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.text.y = element_text(
      size = 11
    ),
    axis.text.x = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    )
  )
print(p_other)


#
library(dplyr)
library(ggplot2)
library(stringr)

# =========================================================
# 1. 整理数据
# =========================================================

go_plot <- go_compare %>%
  dplyr::mutate(
    Group = as.character(Group),
    Description = as.character(Description),
    GeneRatio_num = as.numeric(GeneRatio_num),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust),
    
    # 自动把 GO 名称换行
    Description_wrap = stringr::str_wrap(
      Description,
      width = 30
    )
  )


# =========================================================
# 2. 看看真实的 Group 名称
# =========================================================

print(unique(go_plot$Group))
print(table(go_plot$Group))


# =========================================================
# 3. 统一 Group 显示名称
# =========================================================

go_plot <- go_plot %>%
  dplyr::mutate(
    Group_plot = dplyr::case_when(
      
      grepl(
        "covid",
        Group,
        ignore.case = TRUE
      ) ~ "COVID19-specific",
      
      grepl(
        "other.*virus|virus.*specific",
        Group,
        ignore.case = TRUE
      ) ~ "Other-virus-specific",
      
      TRUE ~ Group
    )
  )


# =========================================================
# 4. 固定两个 panel 的顺序
# =========================================================

go_plot$Group_plot <- factor(
  go_plot$Group_plot,
  levels = c(
    "COVID19-specific",
    "Other-virus-specific"
  )
)


# =========================================================
# 5. 画图
# =========================================================

p <- ggplot(
  go_plot,
  aes(
    x = GeneRatio_num,
    y = reorder(
      Description_wrap,
      GeneRatio_num
    )
  )
) +
  
  geom_point(
    aes(
      size = Count,
      color = p.adjust
    ),
    alpha = 0.85
  ) +
  
  # -------------------------
# 气泡大小
# -------------------------

scale_size_continuous(
  name = "Gene Count",
  range = c(3, 10)
) +
  
  # -------------------------
# p value 颜色
# -------------------------

scale_color_gradient(
  name = "Adjusted p-value",
  low = "#E64B35",
  high = "#4C78A8"
) +
  
  # -------------------------
# X轴
# -------------------------

scale_x_continuous(
  expand = expansion(
    mult = c(0.05, 0.10)
  )
) +
  
  # -------------------------
# 两个 panel
# -------------------------

facet_wrap(
  ~ Group_plot,
  nrow = 1,
  scales = "free_x"
) +
  
  labs(
    title = "GO Biological Process Comparison",
    x = "Gene Ratio",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    
    # 总标题
    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        bottom = 20
      )
    ),
    
    # panel 标题
    strip.text = element_text(
      size = 18,
      face = "bold",
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    # GO term
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    # x轴
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 14,
      margin = margin(
        top = 10
      )
    ),
    
    # 网格
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    
    # panel之间距离
    panel.spacing = unit(
      1.2,
      "cm"
    ),
    
    # legend
    legend.position = "right",
    
    legend.title = element_text(
      size = 14
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

p

##conclusion:COVID19-specific genes were predominantly enriched in antiviral and
##innate immune-related biological processes, whereas genes shared with/associated 
##with other viruses showed stronger enrichment in RNA splicing, autophagy, 
##vesicle transport, and intracellular protein localization processes.

#
df <- read.csv(file.choose(), stringsAsFactors = FALSE)

colnames(df)
dim(df)
table(df$significance, useNA = "ifany")

rm(df)

df <- read.csv(file.choose(), stringsAsFactors = FALSE)

colnames(df)
dim(df)
head(df)

#
dim(df)
colnames(df)
table(df$direction)


df <- read.csv(file.choose(), stringsAsFactors = FALSE)
dim(df)
colnames(df)
table(df$direction)
#

up_genes <- df$gene[df$direction == "Up"]
down_genes <- df$gene[df$direction == "Down"]
length(up_genes)
length(down_genes)
head(up_genes)
head(down_genes)



#drug re

summary(df$log2FoldChange)
summary(df$padj)
head(
  df[order(df$padj), c("gene", "log2FoldChange", "pvalue", "padj", "direction")],
  20
)

#
table(
df$direction,
df$padj < 0.05,
useNA = "ifany"
)
table(
  df$direction,
  abs(df$log2FoldChange) > 1,
  useNA = "ifany"
)
summary(df$log2FoldChange[df$direction == "Up"])
summary(df$log2FoldChange[df$direction == "Down"])



#
write.table(
  up_genes,
  "COVID19_UP_signature.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  down_genes,
  "COVID19_DOWN_signature.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
length(up_genes)
length(down_genes)

write.table(up_genes, "COVID19_UP_signature.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(down_genes, "COVID19_DOWN_signature.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

getwd()



#top100
top_up <- df[df$direction == "Up", ]
top_up <- top_up[order(top_up$padj), ]
top_up <- head(top_up, 100)
top_down <- df[df$direction == "Down", ]
top_down <- top_down[order(top_down$padj), ]
top_down <- head(top_down, 100)
nrow(top_up)
nrow(top_down)
head(top_up[, c("gene", "log2FoldChange", "padj")])
head(top_down[, c("gene", "log2FoldChange", "padj")])

write.table(
  top_up[, "gene", drop = FALSE],
  "COVID19_Top100_UP.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  top_down[, "gene", drop = FALSE],
  "COVID19_Top100_DOWN.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
getwd()


###
file.exists("COVID19_Top100_UP.txt")
file.exists("COVID19_Top100_DOWN.txt")
length(readLines("COVID19_Top100_UP.txt"))
length(readLines("COVID19_Top100_DOWN.txt"))
