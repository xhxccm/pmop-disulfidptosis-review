# ============================================================
# 用 R 原生包重做其余 Python 绘制的图（ggplot2/pheatmap/VennDiagram）
#   Fig3B SVM-RFE 重要性  Fig3C RF 重要性  Fig3D 三算法韦恩
#   Fig5A 免疫堆叠  Fig5B Hub-免疫相关热图  Fig6C 对接结合能柱状图
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap); library(VennDiagram); library(grid)
})
OUT <- "figures/R_native"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- Fig3B SVM-RFE 重要性 ----------
svm <- read.csv("hub_genes_results/SVMRFE_feature_ranking.csv", row.names = 1)
svm$gene <- rownames(svm); svm$ranking <- as.integer(svm$ranking)
svm_top <- head(svm[order(svm$ranking), ], 30)
svm_top$importance <- max(svm_top$ranking) - svm_top$ranking + 1
p <- ggplot(svm_top, aes(x = reorder(gene, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "#4A7FB5") + coord_flip() +
  labs(title = "SVM-RFE feature ranking (top 30)", x = "Gene", y = "Relative importance") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUT, "Fig3B_SVMRFE_importance_R.png"), p, width = 8, height = 7, dpi = 200)

# ---------- Fig3C RF 重要性（Gini + 置换双指标）----------
rf <- read.csv("hub_genes_results/RF_feature_importance.csv", row.names = 1)
rf$gene <- rownames(rf)
rf_top <- head(rf[order(-rf$gini_importance), ], 25)
rf_long <- data.frame(gene = rep(rf_top$gene, 2),
                      metric = rep(c("Gini", "Permutation"), each = nrow(rf_top)),
                      value = c(rf_top$gini_importance, rf_top$perm_importance))
p <- ggplot(rf_long, aes(x = reorder(gene, value), y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") + coord_flip() +
  facet_wrap(~ metric, scales = "free_x") +
  scale_fill_manual(values = c(Gini = "#d62728", Permutation = "#1f77b4")) +
  labs(title = "Random Forest feature importance (top 25)", x = "Gene", y = "Importance") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
ggsave(file.path(OUT, "Fig3C_RF_importance_R.png"), p, width = 9, height = 8, dpi = 200)

# ---------- Fig3D 三算法交集韦恩（用 vote 布尔列，与 Python 版一致）----------
vote <- read.csv("hub_genes_results/three_algorithm_vote.csv", check.names = FALSE)
tobool <- function(x) as.character(x) %in% c("True", "TRUE", "true", "1")
lasso <- vote$gene_symbol[tobool(vote$in_LASSO)]
svm_g <- vote$gene_symbol[tobool(vote$in_SVMRFE)]
rf_g  <- vote$gene_symbol[tobool(vote$in_RF)]
venn.diagram(list(LASSO = lasso, `SVM-RFE` = svm_g, `Random Forest` = rf_g),
             filename = file.path(OUT, "Fig3D_three_algorithm_venn_R.png"),
             fill = c("#4A9E9E", "#4A7FB5", "#B57A4A"), alpha = 0.5, cex = 1.1, cat.cex = 1.0,
             main = "LASSO / SVM-RFE / RF intersection", main.cex = 1.2,
             height = 1800, width = 1800)
cat(sprintf("3-algo intersection: %s\n", paste(Reduce(intersect, list(lasso, svm_g, rf_g)), collapse = ", ")))

# ---------- Fig5A 免疫细胞堆叠图 ----------
frac <- read.csv("immune_results/immune_fraction.csv", row.names = 1, check.names = FALSE)
info <- read.csv("GSE56815_sample_info.csv", row.names = 1, check.names = FALSE)
info <- info[rownames(frac), , drop = FALSE]
ord <- order(info$group); frac <- frac[ord, , drop = FALSE]; info <- info[ord, , drop = FALSE]
celltypes <- colnames(frac)          # 先固定 22 类细胞名，勿把 sample 列算进去
samples <- rownames(frac)
fl <- data.frame(sample = rep(samples, length(celltypes)),
                 celltype = rep(celltypes, each = length(samples)),
                 fraction = as.vector(as.matrix(frac[, celltypes, drop = FALSE])))
fl$group <- info$group[match(fl$sample, rownames(info))]
fl$sample <- factor(fl$sample, levels = samples)
p <- ggplot(fl, aes(x = sample, y = fraction, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +
  labs(title = "Immune cell infiltration (CIBERSORT-style NNLS)", x = "Sample", y = "Fraction", fill = "Cell type") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                     axis.text.x = element_blank(), legend.position = "right",
                     legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))
ggsave(file.path(OUT, "Fig5A_immune_stacked_R.png"), p, width = 12, height = 6, dpi = 200)
cat("[STEP] Fig5A done\n")

# ---------- Fig5B Hub-免疫相关热图 ----------
corr <- as.matrix(read.csv("immune_results/HubGene_immune_correlation.csv", row.names = 1, check.names = FALSE))
cat("[STEP] Fig5B pheatmap start; NA in corr:", sum(is.na(corr)), "\n")
pheatmap(corr, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#1f77b4", "white", "#d62728"))(100),
         breaks = seq(-1, 1, length.out = 101), display_numbers = FALSE,
         fontsize_row = 7, fontsize_col = 6, angle_col = 45,
         main = "Hub gene - immune cell correlation (Spearman)",
         filename = file.path(OUT, "Fig5B_immune_correlation_R.png"), width = 12, height = 8)
cat("Fig3B/3C/3D + Fig5A/5B (R) done\n")

# ---------- Fig6C 对接结合能柱状图 ----------
dock <- read.csv("docking/results/docking_affinity_matrix.csv", row.names = 1, check.names = FALSE)
best <- read.csv("docking/results/best_ligand_per_receptor.csv")
# 取每个受体的最优结合能（矩阵每列最小值）
best_ene <- data.frame(receptor = colnames(dock),
                       affinity = sapply(dock, min),
                       ligand = sapply(colnames(dock), function(cn) rownames(dock)[which.min(dock[[cn]])]))
best_ene$receptor <- sub("_.*", "", best_ene$receptor)
p <- ggplot(best_ene, aes(x = reorder(receptor, affinity), y = affinity, fill = affinity)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_gradient(low = "#d62728", high = "#1f77b4") +
  geom_text(aes(label = sprintf("%.2f (%s)", affinity, ligand)), hjust = -0.05, size = 3) +
  labs(title = "Molecular docking: best binding affinity per receptor",
       x = "Receptor", y = "Binding affinity (kcal/mol)") +
  ylim(min(best_ene$affinity) * 1.05, 0) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
ggsave(file.path(OUT, "Fig6C_docking_bar_R.png"), p, width = 9, height = 6, dpi = 200)
cat("Fig6C (R) done. ALL DONE\n")
