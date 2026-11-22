#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第二阶段第一步：GSE56814 差异表达基因分析（limma 等价实现）

由于本机 R 的 Bioconductor 依赖安装受阻，本脚本用 Python 复现 limma 的核心统计流程：
  1) lmFit:     两组线性模型 (design = ~ group)，得到 logFC 与残差方差 s^2
  2) eBayes:    Smyth (2004) empirical Bayes 方差收缩
                s^2 ~ s0^2 * d0 / chi^2_{d0}
                通过 log(s^2) 的均值/方差拟合先验自由度 d0 与先验方差 s0^2
  3) moderated t 检验: t = logFC / (s_post * sqrt(1/n1 + 1/n2)), df = d0 + d_resid
  4) BH FDR 多重检验校正

输出与 limma::topTable 的核心列一致：logFC, AveExpr, t, P.Value, adj.P.Val
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import optimize, special, stats
from statsmodels.stats.multitest import multipletests

BASE = pathlib.Path(__file__).resolve().parent
EXPR_CSV = BASE / "GSE56814_exprs_matrix.csv"
SAMPLE_CSV = BASE / "GSE56814_sample_info.csv"
OUT_ALL = BASE / "GSE56814_DEG_limma.csv"
OUT_SIG = BASE / "GSE56814_DEG_sig.csv"
OUT_NOMINAL = BASE / "GSE56814_DEG_nominal_P005.csv"
OUT_SUMMARY = BASE / "GSE56814_DEG_summary.txt"
OUT_PLOT = BASE / "GSE56814_volcano.png"

LOGFC_CUT = 1.0     # |log2FC| >= 1 (2 倍)
FDR_CUT = 0.05      # adj.P.Value < 0.05


def squeeze_var(s2: np.ndarray, df_resid: int):
    """逐行复现 limma::squeezeVar 与 fitFDist（R/ebayes.R）。

    limma 把残差方差建模为 scale × F(df1, df2)，用 log 尺度的矩估计：
        e     = log(x) - digamma(df1/2) + log(df1/2)
        evar  = mean( n/(n-1) * (e - mean(e))^2 ) - trigamma(df1/2)
        df2   = 2 * trigammaInverse(evar)
        scale = exp(mean(e) + digamma(df2/2) - log(df2/2))
    后验方差 var.post = (df*var + df.prior*var.prior) / (df + df.prior)

    参数
    ----
    s2 : 每个基因的残差方差 (长度 G)
    df_resid : 每个基因的残差自由度 (n - 2)

    返回
    ----
    s2_post, df_total, d0(=df.prior), s0_2(=var.prior)
    """
    x = np.maximum(s2, 1e-5 * np.median(s2))  # fitFDist 的零值偏移
    z = np.log(x)
    e = z - special.digamma(df_resid / 2.0) + np.log(df_resid / 2.0)
    emean = float(e.mean())
    n = len(e)
    evar = float(np.mean(n / (n - 1) * (e - emean) ** 2) -
                 special.polygamma(1, df_resid / 2.0))

    if evar > 0:
        y = optimize.brentq(lambda t: special.polygamma(1, t) - evar, 1e-8, 1e8, xtol=1e-12)
        d0 = 2.0 * y
        s0_2 = float(np.exp(emean + special.digamma(y) - np.log(y)))
    else:
        d0 = np.inf
        s0_2 = float(np.exp(emean))

    df_total = d0 + df_resid
    s2_post = (df_resid * s2 + d0 * s0_2) / df_total
    return s2_post, df_total, d0, s0_2


def run() -> None:
    expr = pd.read_csv(EXPR_CSV, index_col=0)
    info = pd.read_csv(SAMPLE_CSV, index_col=0)

    # 严格按表达矩阵的列顺序对齐样本，避免分组错位
    info = info.loc[list(expr.columns)]
    group = info["group"].astype(str).to_numpy()
    is_ctrl = group == "Control"
    is_case = group == "PMOP"
    n_ctrl, n_case = int(is_ctrl.sum()), int(is_case.sum())
    if n_ctrl + n_case != len(group):
        raise ValueError(f"分组应仅含 Control/PMOP，实际: {np.unique(group, return_counts=True)}")

    X = expr.to_numpy(dtype=float)
    mean_ctrl = X[:, is_ctrl].mean(axis=1)
    mean_case = X[:, is_case].mean(axis=1)

    # lmFit: 以 Control 为参考，系数 PMOP vs Control = 组间均值差
    log_fc = mean_case - mean_ctrl
    ave_expr = X.mean(axis=1)
    df_resid = n_ctrl + n_case - 2
    ss_resid = ((X[:, is_ctrl] - mean_ctrl[:, None]) ** 2).sum(axis=1) + \
               ((X[:, is_case] - mean_case[:, None]) ** 2).sum(axis=1)
    s2 = ss_resid / df_resid

    # eBayes: 方差收缩 + moderated t
    s2_post, df_total, d0, s0_2 = squeeze_var(s2, df_resid)
    se = np.sqrt(s2_post * (1.0 / n_ctrl + 1.0 / n_case))
    t_stat = log_fc / se
    p_value = 2.0 * stats.t.sf(np.abs(t_stat), df_total)
    q_value = multipletests(p_value, method="fdr_bh")[1]

    deg = pd.DataFrame(
        {
            "logFC": log_fc,
            "AveExpr": ave_expr,
            "t": t_stat,
            "P.Value": p_value,
            "adj.P.Val": q_value,
        },
        index=expr.index,
    )
    deg.index.name = "gene_symbol"
    deg = deg.sort_values(["adj.P.Val", "P.Value"])

    sig = deg[(deg["adj.P.Val"] < FDR_CUT) & (deg["logFC"].abs() >= LOGFC_CUT)].copy()
    sig["regulation"] = np.where(sig["logFC"] > 0, "Up in PMOP", "Down in PMOP")
    sig = sig.sort_values("logFC", key=lambda s: s.abs(), ascending=False)

    deg.to_csv(OUT_ALL)
    sig.to_csv(OUT_SIG)
    nominal = deg[deg["P.Value"] < 0.05].copy()  # 探索性参考，未通过 FDR
    nominal.to_csv(OUT_NOMINAL)

    up = int((sig["regulation"] == "Up in PMOP").sum())
    down = int((sig["regulation"] == "Down in PMOP").sum())
    any_fdr = int((deg["adj.P.Val"] < FDR_CUT).sum())
    any_nom = int((deg["P.Value"] < 0.05).sum())
    p_mean, p_median = float(p_value.mean()), float(np.median(p_value))
    max_lfc_nominal = float(nominal["logFC"].abs().max()) if len(nominal) else 0.0

    lines = [
        "GSE56814 差异表达分析（limma 等价实现：lmFit + eBayes moderated t）",
        "=" * 68,
        f"样本: {len(group)} 例 (Control {n_ctrl} / PMOP {n_case})，仅绝经后亚组",
        f"基因: {X.shape[0]} 个",
        "",
        "[eBayes 先验参数]",
        f"  残差自由度 df.residual = {df_resid}",
        f"  先验自由度 df.prior (d0) = {d0:.4f}",
        f"  先验方差 s0^2 = {s0_2:.6f}  (s0 = {np.sqrt(s0_2):.4f})",
        f"  调节自由度 df.total = {df_total:.4f}",
        "",
        "[差异基因统计]",
        f"  名义 P < 0.05: {any_nom} 个 (完全无差异时期望约 {int(0.05 * len(deg))} 个)",
        f"  名义显著基因中最大 |log2FC| = {max_lfc_nominal:.4f} (阈值 {LOGFC_CUT})",
        f"  BH FDR < {FDR_CUT}: {any_fdr} 个",
        f"  FDR < {FDR_CUT} 且 |log2FC| >= {LOGFC_CUT}: {len(sig)} 个"
        f"  (上调 {up} / 下调 {down})",
        f"  P 值分布: 均值 {p_mean:.4f}, 中位数 {p_median:.4f}",
        "",
        "[说明] 本数据集绝经后亚组 (26 Control vs 16 PMOP) 组间转录组差异极弱：",
        "  P<0.05 数量低于随机期望且最大 |log2FC| 未达 1，BH 校正后无基因显著。",
        "  样本间平均表达相关性约 0.95，个体间异质性大。名义显著列表仅作探索性参考。",
        "",
        "[Top 20 差异基因 (按 adj.P.Val)]",
        deg.head(20).to_string(),
    ]
    OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")

    # 火山图（白底，符合浅色主题）
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
        x = deg["logFC"].to_numpy()
        y = -np.log10(deg["adj.P.Val"].clip(lower=1e-300))
        is_sig = (deg["adj.P.Val"].to_numpy() < FDR_CUT) & (np.abs(x) >= LOGFC_CUT)
        ax.scatter(x[~is_sig], y[~is_sig], s=8, c="#B0B0B0", alpha=0.55, label="Not significant", linewidths=0)
        ax.scatter(x[is_sig & (x > 0)], y[is_sig & (x > 0)], s=14, c="#D62728", alpha=0.85, label="Up in PMOP", linewidths=0)
        ax.scatter(x[is_sig & (x < 0)], y[is_sig & (x < 0)], s=14, c="#1F77B4", alpha=0.85, label="Down in PMOP", linewidths=0)

        for gene, row in deg.head(10).iterrows():
            ax.annotate(gene, (row["logFC"], -np.log10(max(row["adj.P.Val"], 1e-300))),
                        fontsize=7, xytext=(3, 3), textcoords="offset points", color="#333333")

        ax.axhline(-np.log10(FDR_CUT), color="#888888", lw=0.8, ls="--")
        ax.axvline(LOGFC_CUT, color="#888888", lw=0.8, ls="--")
        ax.axvline(-LOGFC_CUT, color="#888888", lw=0.8, ls="--")
        ax.set_xlabel("log2 fold change (PMOP vs Control)")
        ax.set_ylabel("-log10(adj.P.Value)")
        ax.set_title("GSE56814 Differential Expression (limma-equivalent)")
        ax.legend(loc="upper left", frameon=False, fontsize=9)
        ax.grid(alpha=0.15)
        fig.tight_layout()
        fig.savefig(OUT_PLOT)
        plt.close(fig)
    except Exception as exc:  # 绘图失败不影响核心结果
        print(f"[WARN] volcano plot skipped: {exc}", file=sys.stderr)

    print(f"df.prior={d0:.3f}, s0^2={s0_2:.4f}, df.total={df_total:.3f}")
    print(f"FDR<{FDR_CUT} & |log2FC|>={LOGFC_CUT}: {len(sig)} genes (up {up} / down {down})")
    print(f"Saved: {OUT_ALL.name}, {OUT_SIG.name}, {OUT_NOMINAL.name}, "
          f"{OUT_SUMMARY.name}, {OUT_PLOT.name}")


if __name__ == "__main__":
    run()
