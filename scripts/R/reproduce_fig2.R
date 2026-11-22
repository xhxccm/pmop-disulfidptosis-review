# ============================================================
# 用 R 原生包重做 Fig.2 系列（此前因 R 不可用而用 Python 绘制）
#   Fig2A 火山图      -> ggplot2
#   Fig2B DEG 热图    -> pheatmap
#   Fig2C DRG-DEG 韦恩 -> VennDiagram
# 数据全部来自既有分析产物，不重新计算统计
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap); library(VennDiagram); library(grid)
})
# 工作目录由调用方（Bash）继承，勿在脚本内 setwd 中文路径（C locale 下会乱码）
OUT <- "figures/R_native"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

deg   <- read.csv("GSE56815_DEG_limma.csv", stringsAsFactors = FALSE)
expr  <- read.csv("GSE56815_exprs_matrix.csv", row.names = 1, check.names = FALSE)
info  <- read.csv("GSE56815_sample_info.csv", row.names = 1, check.names = FALSE)
drg   <- read.csv("disulfidptosis_geneset_extended.csv")$gene_symbol

# ---------- Fig2A 火山图 ----------
deg$sig <- "ns"
deg$sig[deg$P.Value < 0.05 & deg$adj.P.Val >= 0.05 & deg$logFC > 0] <- "Nominal_up"
deg$sig[deg$P.Value < 0.05 & deg$adj.P.Val >= 0.05 & deg$logFC < 0] <- "Nominal_dn"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC > 0] <- "FDR_up"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC < 0] <- "FDR_dn"
deg$sig <- factor(deg$sig, levels = c("ns", "Nominal_up", "Nominal_dn", "FDR_up", "FDR_dn"))
cols <- c(ns = "grey80", Nominal_up = "#f6b26b", Nominal_dn = "#9fc5e8",
          FDR_up = "#d62728", FDR_dn = "#1f77b4")
p <- ggplot(deg, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(size = 1.1, alpha = 0.8) +
  scale_color_manual(values = cols, labels = c("ns", "Nominal up", "Nominal down", "FDR up", "FDR down")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "GSE56815: Differential expression (limma)",
       x = "log2 Fold Change", y = "-log10(P.Value)", color = "Group") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                     legend.position = "top")
# 标注 Top10 最显著基因
top10 <- deg[order(deg$adj.P.Val)[1:10], ]
p <- p + geom_text(data = top10, aes(label = gene_symbol), size = 2.8, color = "black", vjust = -0.6)
ggsave(file.path(OUT, "Fig2A_volcano_R.png"), p, width = 8, height = 6, dpi = 200)
cat("Fig2A volcano (R) done\n")

# ---------- Fig2B DEG 热图 ----------
sig <- deg[order(deg$adj.P.Val), ]
sig <- sig$gene_symbol[1:60]
sig <- intersect(sig, rownames(expr))
mat <- as.matrix(expr[sig, , drop = FALSE])
# 每行 Z-score
mat <- t(scale(t(mat)))
mat[mat > 2.5] <- 2.5; mat[mat < -2.5] <- -2.5
# 按 group 排序样本
info <- info[colnames(mat), , drop = FALSE]
ord <- order(info$group)
mat <- mat[, ord, drop = FALSE]
ann <- data.frame(Group = info$group[ord], row.names = colnames(mat))
ann_col <- list(Group = c(Control = "#1f77b4", PMOP = "#d62728"))
pheatmap(mat, annotation_col = ann, annotation_colors = ann_col,
         show_rownames = TRUE, show_colnames = FALSE, cluster_cols = FALSE,
         cluster_rows = TRUE, fontsize_row = 6,
         color = colorRampPalette(c("#1f77b4", "white", "#d62728"))(100),
         main = "Top 60 FDR-significant DEGs (Z-score)",
         filename = file.path(OUT, "Fig2B_DEG_heatmap_R.png"),
         width = 7, height = 8)
cat("Fig2B heatmap (R) done\n")

# ---------- Fig2C DRG-DEG 韦恩图 ----------
deg_fdr  <- deg$gene_symbol[deg$adj.P.Val < 0.05]
deg_nom  <- deg$gene_symbol[deg$P.Value < 0.05]
venn.diagram(list(A = drg, B = deg_fdr), filename = file.path(OUT, "Fig2C_DRG_FDRdeg_venn_R.png"),
             fill = c("#5B9E78", "#4A7FB5"), alpha = 0.5, cex = 1.1, cat.cex = 1.0,
             category.names = c(paste0("DRGs\n(n=", length(drg), ")"),
                                paste0("FDR DEGs\n(n=", length(deg_fdr), ")")),
             main = "DRGs vs FDR-significant DEGs", main.cex = 1.2, height = 1500, width = 1500)
venn.diagram(list(A = drg, B = deg_nom), filename = file.path(OUT, "Fig2C_DRG_Nominaldeg_venn_R.png"),
             fill = c("#5B9E78", "#B57A4A"), alpha = 0.5, cex = 1.1, cat.cex = 1.0,
             category.names = c(paste0("DRGs\n(n=", length(drg), ")"),
                                paste0("Nominal DEGs\n(n=", length(deg_nom), ")")),
             main = "DRGs vs Nominal DEGs (P<0.05)", main.cex = 1.2, height = 1500, width = 1500)
cat("Fig2C Venn (R) done\n")
cat(sprintf("Intersection FDR: %s\n", paste(intersect(drg, deg_fdr), collapse = ", ")))
cat(sprintf("Intersection Nominal: %s\n", paste(intersect(drg, deg_nom), collapse = ", ")))
cat("ALL DONE\n")
