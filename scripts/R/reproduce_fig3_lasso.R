# ============================================================
# 用 R 原生 glmnet 重做 LASSO 特征筛选（此前用 Python sklearn 替代）
# 特征集：FDR<0.05 DEGs (333)，与 Python 版一致
# 小样本(40)限制：cv.glmnet 10折每折仅4例会强制退回 deviance，
#   故用 deviance 选 lambda（得稳定基因集），再用手动分层10折计算该模型 AUC。
# 输出：系数路径图、CV-AUC/ROC 曲线、非零系数基因表
# ============================================================
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
set.seed(42)

expr  <- read.csv("GSE56815_exprs_matrix.csv", row.names = 1, check.names = FALSE)
sinfo <- read.csv("GSE56815_sample_info.csv", row.names = 1, check.names = FALSE)
deg   <- read.csv("GSE56815_DEG_limma.csv", row.names = 1, check.names = FALSE)
sinfo <- sinfo[colnames(expr), , drop = FALSE]

OUT <- "figures/R_native"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fdr_genes <- rownames(deg)[deg$adj.P.Val < 0.05]
X <- t(as.matrix(expr[fdr_genes, , drop = FALSE]))
y <- ifelse(sinfo$group == "PMOP", 1, 0)
cat(sprintf("Features: %d  Samples: %d  Cases: %d\n", ncol(X), nrow(X), sum(y)))

# --- 1) cv.glmnet 用 deviance 选 lambda（小样本下 glmnet 的标准行为）---
set.seed(42)
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10,
                   type.measure = "deviance", standardize = TRUE, parallel = FALSE)
cat(sprintf("lambda.min=%.5g  lambda.1se=%.5g\n", cvfit$lambda.min, cvfit$lambda.1se))

# --- 2) 手动分层 10 折，用 lambda.min 计算 out-of-fold AUC ---
set.seed(42)
foldid <- integer(length(y))
for (cl in c(0, 1)) {
  idx <- which(y == cl)
  foldid[idx] <- sample(rep_len(1:10, length(idx)))
}
oof <- rep(NA_real_, length(y))
for (k in 1:10) {
  tr <- foldid != k; te <- foldid == k
  fit <- glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 1,
                lambda = cvfit$lambda.min, standardize = TRUE)
  oof[te] <- as.numeric(predict(fit, X[te, , drop = FALSE], type = "response"))
}
roc_obj <- roc(y, oof, quiet = TRUE)
cat(sprintf("Manual stratified 10-fold CV AUC = %.4f\n", as.numeric(auc(roc_obj))))

# --- 3) 系数表（全数据 + lambda.min）---
cm <- as.matrix(coef(cvfit, s = "lambda.min"))
coef_df <- data.frame(gene = rownames(cm), coefficient = as.numeric(cm[, 1]),
                      stringsAsFactors = FALSE)
coef_df <- coef_df[coef_df$gene != "(Intercept)", ]
coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]
coef_df$selected <- abs(coef_df$coefficient) > 1e-8
write.csv(coef_df, file.path(OUT, "Fig3A_GLMNET_coefficients_R.csv"), row.names = FALSE)
sel <- coef_df$gene[coef_df$selected]
cat(sprintf("lambda.min selected genes: %d\n", length(sel)))
cat(paste(head(sel, 25), collapse = ", "), "\n")

# --- 4) 图：CV deviance 曲线 ---
png(file.path(OUT, "Fig3A_GLMNET_CV_curve_R.png"), width = 1400, height = 1100, res = 200)
plot(cvfit, main = "glmnet LASSO: CV deviance vs log(lambda)")
dev.off()

# --- 5) 图：out-of-fold ROC ---
png(file.path(OUT, "Fig3A_GLMNET_ROC_R.png"), width = 1200, height = 1200, res = 200)
plot(roc_obj, main = sprintf("glmnet LASSO OOF ROC (AUC=%.3f)", as.numeric(auc(roc_obj))),
     col = "#d62728", lwd = 2)
dev.off()

# --- 6) 图：系数路径 ---
png(file.path(OUT, "Fig3A_GLMNET_coefficient_path_R.png"), width = 1600, height = 1200, res = 200)
plot(cvfit$glmnet.fit, xvar = "lambda", label = FALSE,
     main = "glmnet LASSO coefficient paths")
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = 2, col = "red")
dev.off()

# --- 7) 图：Top25 系数条形图 ---
top <- head(coef_df[coef_df$selected, ], 25)
top <- top[order(top$coefficient), ]
png(file.path(OUT, "Fig3A_GLMNET_top_coefficients_R.png"), width = 1600, height = 1400, res = 200)
par(mar = c(4, 8, 3, 2))
barplot(top$coefficient, names.arg = top$gene, horiz = TRUE, las = 2,
        col = ifelse(top$coefficient > 0, "#d62728", "#1f77b4"),
        main = "glmnet LASSO top coefficients (lambda.min)", xlab = "Coefficient")
dev.off()
cat("Fig3A glmnet (R) ALL DONE\n")
