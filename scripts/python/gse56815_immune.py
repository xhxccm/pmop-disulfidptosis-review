#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第五阶段 5.1-5.2：免疫细胞浸润（CIBERSORT 近似）+ Hub Gene-免疫相关性

5.1 CIBERSORT 反卷积：内置覆盖 LM22 全部 22 类免疫细胞的经典 marker 签名，
    用 Nu-SVR 做反卷积，输出 22 类细胞比例 stacked bar。
    注意：完整 LM22 为 547 基因签名矩阵，本实现使用各细胞的经典 marker 基因集，
    属近似实现，结果用于趋势比较而非绝对定量。
5.2 Hub Genes 与免疫细胞比例的 Spearman 相关 + 热图。
"""
from __future__ import annotations

import sys
import pathlib
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "py_lib"))
from sklearn.svm import NuSVR
from sklearn.preprocessing import StandardScaler

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "immune_results"
OUT.mkdir(exist_ok=True)
EXPR = BASE / "GSE56815_exprs_matrix.csv"
SAMPLE = BASE / "GSE56815_sample_info.csv"

# LM22 22 类免疫细胞的经典 marker 基因（Newman et al. 2015 及常用免疫图谱）
LM22_MARKERS = {
    "B cells naive": ["TCL1A", "IGHD", "IL4R", "CD79A", "MS4A1", "CD19", "FCER2"],
    "B cells memory": ["CD27", "CD80", "TNFRSF13B", "BANK1", "AIM2", "CD40", "CD79B"],
    "Plasma cells": ["MZB1", "XBP1", "IRF4", "CD38", "PRDM1", "JCHAIN", "SDC1"],
    "T cells CD8": ["CD8A", "CD8B", "GZMK", "CCL5", "GZMB", "NKG7", "CD3D"],
    "T cells CD4 naive": ["CCR7", "TCF7", "LEF1", "SELL", "CD27", "MAL", "IL7R"],
    "T cells CD4 memory resting": ["CD4", "IL7R", "KLF2", "LTB", "CD27", "TXK"],
    "T cells CD4 memory activated": ["HLA-DRA", "CD40LG", "ICOS", "TNFRSF4", "CD69", "CD4"],
    "T cells follicular helper": ["CXCR5", "BCL6", "IL21", "PDCD1", "ICOS", "CXCL13"],
    "T cells regulatory (Tregs)": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "TIGIT"],
    "T cells gamma delta": ["TRDC", "TRGC1", "TRGC2", "CD3D", "CD3E"],
    "NK cells resting": ["KLRD1", "NCAM1", "FCGR3A", "NKG7", "GNLY", "KLRB1"],
    "NK cells activated": ["FCGR3A", "NKG7", "GZMB", "PRF1", "IFNG", "KLRD1"],
    "Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "CSF1R", "FCGR1A", "VCAN"],
    "Macrophages M0": ["CD68", "CD163", "CSF1R", "MRC1", "MSR1", "CD14"],
    "Macrophages M1": ["NOS2", "CD86", "IL1B", "TNF", "CD80", "HLA-DRA"],
    "Macrophages M2": ["MRC1", "CD163", "MS4A7", "ALOX15", "CD209", "STAB1"],
    "Dendritic cells resting": ["CD1C", "FCER1A", "CST3", "CLEC10A", "ITGAX"],
    "Dendritic cells activated": ["CD83", "CLEC9A", "HLA-DRA", "CD40", "IL12B", "CCL22"],
    "Mast cells resting": ["TPSAB1", "CPA3", "MS4A2", "KIT", "HDC"],
    "Mast cells activated": ["TPSAB1", "CPA3", "IL1B", "TNF", "CCL3", "MS4A2"],
    "Eosinophils": ["EPX", "RNASE2", "RNASE3", "CLC", "SIGLEC8", "IL5RA"],
    "Neutrophils": ["FCGR3B", "CSF3R", "CXCR1", "CXCR2", "MPO", "ELANE", "S100A8"],
}


def deconvolve(expr, markers):
    """CIBERSORT 式线性反卷积：V(基因x样本) ≈ B(基因x细胞类型) @ P(细胞类型x样本)。
    用非负最小二乘(NNLS)逐样本求解 P，得到各细胞类型相对丰度。"""
    from scipy.optimize import nnls
    genes = sorted({g for gl in markers.values() for g in gl} & set(expr.index))
    ctypes = list(markers.keys())
    n_types = len(ctypes)
    gene_idx = {g: i for i, g in enumerate(genes)}
    # 参考签名矩阵 B: 基因 x 细胞类型，marker=1，其余=0，按列归一化
    B = np.zeros((len(genes), n_types))
    for t, (ctype, gl) in enumerate(markers.items()):
        for g in gl:
            if g in gene_idx:
                B[gene_idx[g], t] = 1.0
    B = B / (np.linalg.norm(B, axis=0, keepdims=True) + 1e-8)
    # 表达矩阵按基因 z-score（跨样本），与签名可比
    V = expr.loc[genes].values
    V = (V - V.mean(axis=1, keepdims=True)) / (V.std(axis=1, keepdims=True) + 1e-8)
    V = np.clip(V, 0, None)  # NNLS 要求非负
    fractions = np.zeros((V.shape[1], n_types))
    for s in range(V.shape[1]):
        P, _ = nnls(B, V[:, s])
        fractions[s, :] = P
    frac_df = pd.DataFrame(fractions, columns=ctypes, index=expr.columns)
    row_sum = frac_df.sum(axis=1).replace(0, 1)
    frac_df = frac_df.div(row_sum, axis=0)
    return frac_df, genes


def plot_stacked(frac_df, sample, fname):
    order = frac_df.mean().sort_values(ascending=False).index
    f = frac_df[order]
    fig, ax = plt.subplots(figsize=(11, 4.5))
    x = np.arange(f.shape[0])
    bottom = np.zeros(f.shape[0])
    cmap = plt.cm.tab20(np.linspace(0, 1, f.shape[1]))
    for i, col in enumerate(f.columns):
        ax.bar(x, f[col], bottom=bottom, width=1.0, color=cmap[i], label=col, linewidth=0)
        bottom += f[col].values
    groups = sample.loc[f.index, "group"].values
    ax.set_xticks([])
    ax.set_xlabel(f"Samples (first {int((groups=='Control').sum())} Control, rest PMOP)")
    ax.set_ylabel("Estimated fraction")
    ax.set_title("Immune cell infiltration (CIBERSORT-like deconvolution)")
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1), fontsize=7, ncol=1)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=130, bbox_inches="tight")
    plt.close(fig)


def hub_immune_corr(expr, sample, frac_df):
    hubf = BASE / "hub_genes_results" / "HubGenes_three_algorithm_intersection.csv"
    if hubf.exists():
        genes = pd.read_csv(hubf)["gene_symbol"].tolist()
    else:
        genes = pd.read_csv(BASE / "hub_genes_results" / "three_algorithm_vote.csv").head(8)["gene_symbol"].tolist()
    genes = [g for g in genes if g in expr.index]
    if not genes:
        return None, []
    X = expr.loc[genes].T
    corr = pd.DataFrame(index=genes, columns=frac_df.columns, dtype=float)
    pval = pd.DataFrame(index=genes, columns=frac_df.columns, dtype=float)
    for g in genes:
        for c in frac_df.columns:
            r, p = spearmanr(X[g].values, frac_df[c].values)
            corr.loc[g, c] = r
            pval.loc[g, c] = p
    corr.to_csv(OUT / "HubGene_immune_correlation.csv")
    pval.to_csv(OUT / "HubGene_immune_corr_pvalue.csv")
    # 热图
    fig, ax = plt.subplots(figsize=(10, max(3, 0.45 * len(genes))))
    data = corr.values.astype(float)
    im = ax.imshow(data, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data[i,j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("Hub Gene - immune cell Spearman correlation")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(OUT / "HubGene_immune_heatmap.png", dpi=130)
    plt.close(fig)
    return corr, genes


def main():
    expr = pd.read_csv(EXPR, index_col=0)
    sample = pd.read_csv(SAMPLE, index_col=0)
    print(f"表达矩阵: {expr.shape[0]} genes x {expr.shape[1]} samples")
    frac_df, genes = deconvolve(expr, LM22_MARKERS)
    frac_df.to_csv(OUT / "immune_fraction.csv")
    print(f"反卷积完成: {frac_df.shape[1]} 类细胞, marker 基因 {len(genes)} 个")
    print("平均比例(前8):")
    print(frac_df.mean().sort_values(ascending=False).head(8).to_string())
    plot_stacked(frac_df, sample, "immune_infiltration_stacked.png")

    # 组间差异（Wilcoxon 近似用 t 检验）
    ctrl = frac_df[sample.loc[frac_df.index, "group"] == "Control"]
    pmop = frac_df[sample.loc[frac_df.index, "group"] == "PMOP"]
    from scipy.stats import ttest_ind
    diff = pd.DataFrame({
        "control_mean": ctrl.mean(), "pmop_mean": pmop.mean(),
    })
    diff["log2FC"] = np.log2((diff["pmop_mean"] + 1e-6) / (diff["control_mean"] + 1e-6))
    diff["pvalue"] = [ttest_ind(pmop[c], ctrl[c]).pvalue for c in diff.index]
    diff = diff.sort_values("pvalue")
    diff.to_csv(OUT / "immune_fraction_differential.csv")
    print("\n组间差异最显著的前6类:")
    print(diff.head(6).to_string())

    corr, hg = hub_immune_corr(expr, sample, frac_df)
    print(f"\nHub-免疫相关性完成: {len(hg)} genes")
    print("完成。")


if __name__ == "__main__":
    main()
