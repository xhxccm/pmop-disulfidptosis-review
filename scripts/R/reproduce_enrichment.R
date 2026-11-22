# ============================================================
# 用 R 原生 ggplot2 重做富集分析可视化（此前用 Python matplotlib 绘制）
# ORA 气泡图 + GSEA NES 条形图，数据来自 enrichment_results_GSE56815/
# 注：统计值来自既有 ORA/GSEA 计算（Python 版，方法等价），
#     此处仅把绘图引擎换为 R 原生 ggplot2。
# ============================================================
suppressPackageStartupMessages({ library(ggplot2) })
OUT <- "figures/R_native"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RES <- "enrichment_results_GSE56815"

short_term <- function(x) {
  x <- gsub(" \\(GO:[0-9]+\\)", "", x)
  ifelse(nchar(x) > 42, paste0(substr(x, 1, 40), "..."), x)
}

# ---------- ORA 气泡图（扩展 DRG，GO/KEGG）----------
for (lib in c("GO_BP", "GO_CC", "GO_MF", "KEGG")) {
  f <- file.path(RES, sprintf("ORA_DRG_EXTENDED_53_%s.csv", lib))
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- head(d[order(d$adj.P.Val), ], 15)
  d$term <- factor(short_term(d$term), levels = rev(short_term(d$term)))
  p <- ggplot(d, aes(x = fold_enrichment, y = term, size = overlap, color = -log10(adj.P.Val))) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#f6b26b", high = "#d62728", name = expression(-log10(adj.P))) +
    scale_size(name = "Overlap", range = c(3, 10)) +
    labs(title = sprintf("ORA: DRG (extended 53) - %s", lib),
         x = "Fold enrichment", y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                       axis.text.y = element_text(size = 8), legend.position = "right")
  ggsave(file.path(OUT, sprintf("Enrichment_ORA_EXTENDED_%s_R.png", lib)), p, width = 9, height = 6.5, dpi = 200)
}

# ---------- GSEA NES 条形图（KEGG + GO BP）----------
for (lib in c("KEGG", "GO_BP", "GO_CC", "GO_MF")) {
  f <- file.path(RES, sprintf("GSEA_%s.csv", lib))
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) == 0) next
  d <- d[order(d$NES), ]
  d <- head(d, 12)
  d <- rbind(d, head(d[order(-d$NES), ], 12))
  d <- unique(d[, c("term", "NES", "FDR")])
  d <- d[order(d$NES), ]
  d$term <- factor(short_term(d$term), levels = d$term)
  p <- ggplot(d, aes(x = NES, y = term, fill = NES)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#1f77b4", mid = "white", high = "#d62728", midpoint = 0,
                         name = "NES") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = sprintf("GSEA: %s (NES, top up/down)", lib), x = "Normalized ES", y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                       axis.text.y = element_text(size = 7), legend.position = "none")
  ggsave(file.path(OUT, sprintf("Enrichment_GSEA_%s_R.png", lib)), p, width = 9, height = 8, dpi = 200)
}
cat("Enrichment visualization (R ggplot2) ALL DONE\n")
