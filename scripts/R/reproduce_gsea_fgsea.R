# ============================================================
# 用 R 原生 fgsea 重跑 GSEA（此前用 Python prerank 实现）
# 排序指标：limma t 统计量（13038 基因），基因集：本地 GO/KEGG GMT
# 输出 fgsea NES/FDR 表 + 富集图，与 Python GSEA 结果对照
# ============================================================
suppressPackageStartupMessages({ library(fgsea); library(ggplot2) })
set.seed(42)
OUT <- "figures/R_native"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

deg <- read.csv("GSE56815_DEG_limma.csv", stringsAsFactors = FALSE)
stats <- setNames(deg$t, deg$gene_symbol)
stats <- stats[!is.na(stats) & !duplicated(names(stats))]
cat(sprintf("Ranked genes: %d\n", length(stats)))

gmt_files <- c(KEGG = "genesets/KEGG.gmt", GO_BP = "genesets/GO_BP.gmt",
               GO_CC = "genesets/GO_CC.gmt", GO_MF = "genesets/GO_MF.gmt")

for (lib in names(gmt_files)) {
  gl <- gmtPathways(gmt_files[[lib]])
  gl <- gl[sapply(gl, length) >= 15 & sapply(gl, length) <= 500]
  res <- fgseaMultilevel(pathways = gl, stats = stats, minSize = 15, maxSize = 500,
                         nPermSimple = 10000, eps = 1e-10)
  res <- res[order(res$pval), ]
  # leadingEdge 是 list 列，转成字符串才能写 CSV
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
  res$FDR <- res$padj
  write.csv(res, file.path(OUT, sprintf("GSEA_fgsea_%s_R.csv", lib)), row.names = FALSE)
  sig <- res[res$FDR < 0.25, ]
  cat(sprintf("[%s] pathways=%d  FDR<0.25=%d  top: %s\n", lib, nrow(res), nrow(sig),
              paste(head(sig$pathway, 3), collapse = " | ")))

  # 富集分数图：top up + top down
  sel <- c(head(res$pathway[res$NES > 0], 1), head(res$pathway[res$NES < 0], 1))
  sel <- sel[!is.na(sel)]
  if (length(sel) > 0) {
    png(file.path(OUT, sprintf("GSEA_fgsea_%s_enrichment_R.png", lib)),
        width = 1600, height = 1200, res = 200)
    plotEnrichment(gl[[sel[1]]], stats) + labs(title = sprintf("fgsea: %s - %s", lib, sel[1]))
    if (length(sel) > 1) {
      p2 <- plotEnrichment(gl[[sel[2]]], stats) + labs(title = sprintf("fgsea: %s - %s", lib, sel[2]))
      print(gridExtra::grid.arrange(ggplot2::last_plot(), p2, ncol = 1))
    } else {
      print(ggplot2::last_plot())
    }
    dev.off()
  }
}
cat("GSEA fgsea (R) ALL DONE\n")
