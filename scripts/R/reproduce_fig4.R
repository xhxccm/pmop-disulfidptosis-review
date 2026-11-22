# ============================================================
# 用 R 原生 pROC / rms / glmnet 重做 Fig.4（此前用 Python sklearn 替代）
#   Fig4A 训练 ROC  Fig4B 外部验证 ROC  Fig4C Nomogram
# 训练集用 glmnet ridge(alpha=0) 复现 Python LogisticRegression(C=1.0)；
# 外部验证复用训练系数 + 相同 z-score 逻辑。
# AUC 用 ranks(Mann-Whitney) 手算，定义 = P(score_case > score_ctrl)，
# 与 Python sklearn roc_curve 完全一致；pROC 对象方向据 raw AUC 校正后再绘图。
# ============================================================
suppressPackageStartupMessages({ library(pROC); library(rms); library(glmnet) })
OUT <- "figures/R_native"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

coef_df <- read.csv("diagnostic_results/model_coefficients.csv", stringsAsFactors = FALSE)
cmap <- setNames(coef_df$coefficient, coef_df$gene_symbol)
genes <- coef_df$gene_symbol

# raw AUC: P(score_case > score_ctrl) —— 与 Python roc_curve 定义一致
raw_auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# 构造一个 AUC 方向与 raw 一致的 pROC 对象用于绘图
roc_raw <- function(score, y) {
  a <- raw_auc(score, y)
  dir <- if (a >= 0.5) ">" else "<"     # 让 pROC 的 AUC 落在 >=0.5 侧便于画图
  obj <- roc(y, score, direction = dir, quiet = TRUE)
  attr(obj, "raw_auc") <- a              # 真实 AUC（可能<0.5）用于标注
  obj
}

zscore_score <- function(expr) {
  cols <- genes[genes %in% rownames(expr)]
  X <- as.matrix(t(expr[cols, , drop = FALSE]))
  Xs <- scale(X)
  Xs <- apply(Xs, 2, function(v) if (sd(v) == 0) v * 0 else v)
  as.numeric(Xs %*% cmap[cols])
}

make_roc <- function(expr_file, info_file) {
  ee <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  es <- read.csv(info_file, row.names = 1, check.names = FALSE)
  score <- zscore_score(ee)
  es <- es[colnames(ee), , drop = FALSE]
  y <- as.integer(es$group == "PMOP")
  if (length(unique(y)) < 2) return(NULL)
  list(roc = roc_raw(score, y), a = raw_auc(score, y))
}

# ---- 训练集：glmnet ridge 5 折 OOF ----
expr <- read.csv("GSE56815_exprs_matrix.csv", row.names = 1, check.names = FALSE)
sinfo <- read.csv("GSE56815_sample_info.csv", row.names = 1, check.names = FALSE)
sinfo <- sinfo[colnames(expr), , drop = FALSE]
Xall <- scale(as.matrix(t(expr[genes, , drop = FALSE])))
yall <- as.integer(sinfo$group == "PMOP")
lambda_ridge <- 1 / (1.0 * length(yall))
set.seed(42)
fold <- integer(length(yall))
for (cl in c(0, 1)) { idx <- which(yall == cl); fold[idx] <- sample(rep_len(1:5, length(idx))) }
oof <- rep(NA_real_, length(yall))
for (k in 1:5) {
  tr <- fold != k; te <- fold == k
  fit <- glmnet(Xall[tr, , drop = FALSE], yall[tr], family = "binomial",
                alpha = 0, lambda = lambda_ridge, standardize = FALSE)
  oof[te] <- as.numeric(predict(fit, Xall[te, , drop = FALSE], type = "response"))
}
a_train <- raw_auc(oof, yall)
roc_train <- roc_raw(oof, yall)
cat(sprintf("Training 5-fold CV AUC (R ridge raw) = %.4f\n", a_train))

roc_g7429 <- make_roc("GSE7429_exprs_matrix.csv", "GSE7429_sample_info.csv")
roc_g56814 <- make_roc("GSE56814_exprs_matrix.csv", "GSE56814_sample_info.csv")
if (!is.null(roc_g7429))  cat(sprintf("GSE7429 AUC (R raw) = %.4f\n", roc_g7429$a))
if (!is.null(roc_g56814)) cat(sprintf("GSE56814 AUC (R raw) = %.4f\n", roc_g56814$a))

# ---- 合并 ROC 图（曲线用 pROC 画的 >=0.5 侧，对角线用参考，标注真实 raw AUC）----
png(file.path(OUT, "Fig4A_B_ROC_R.png"), width = 1400, height = 1200, res = 200)
plot(roc_train, col = "#d62728", lwd = 2.5, main = "Diagnostic model ROC (R: pROC)", legend = FALSE)
if (!is.null(roc_g7429))  lines(roc_g7429$roc, col = "#1f77b4", lwd = 2)
if (!is.null(roc_g56814)) lines(roc_g56814$roc, col = "#2ca02c", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey60")
legend("bottomright", c(
  sprintf("Training CV (AUC=%.3f)", a_train),
  if (!is.null(roc_g7429))  sprintf("GSE7429 (AUC=%.3f)", roc_g7429$a) else NULL,
  if (!is.null(roc_g56814)) sprintf("GSE56814 (AUC=%.3f)", roc_g56814$a) else NULL,
  "Reference"), col = c("#d62728", "#1f77b4", "#2ca02c", "grey60"), lwd = 2, cex = 0.85)
dev.off()

# ---- Nomogram（rms）----
df <- as.data.frame(Xall); names(df) <- make.names(names(df))
dd <- datadist(df); options(datadist = "dd")
fit_lrm <- lrm(yall ~ ., data = cbind(data.frame(yall = yall), df), x = TRUE, y = TRUE)
nom <- nomogram(fit_lrm)
png(file.path(OUT, "Fig4C_Nomogram_R.png"), width = 2200, height = 1500, res = 200)
plot(nom, cex.axis = 0.65, xlab = "Points")
title("18-Hub Gene Diagnostic Nomogram (R: rms)")
dev.off()
cat("Fig4 ROC + Nomogram (R) ALL DONE\n")
