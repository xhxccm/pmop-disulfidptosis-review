#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第四阶段：诊断模型构建与验证（GSE56815 训练集 + 外部验证集）

4.1 ROC 曲线：训练集 5-fold CV 预测 AUC + 外部集验证 AUC
4.2 列线图 Nomogram：基于训练集逻辑回归系数
4.3 校准曲线（Calibration, 5-fold）+ 决策曲线分析（DCA）
4.4 外部数据集验证：GSE7429（B细胞，同平台）；GSE56814（跨平台，同队列）

Hub Genes 来源：hub_genes_results/HubGenes_three_algorithm_intersection.csv
若交集为空，回退到三算法投票 top N 基因。
"""
from __future__ import annotations

import sys
import pathlib
import warnings
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "py_lib"))
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve, auc, roc_auc_score, brier_score_loss

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "diagnostic_results"
OUT.mkdir(exist_ok=True)
CV = StratifiedKFold(5, shuffle=True, random_state=42)


def load_hub_genes(max_n=10):
    hubf = BASE / "hub_genes_results" / "HubGenes_three_algorithm_intersection.csv"
    if hubf.exists():
        genes = pd.read_csv(hubf)["gene_symbol"].tolist()
        if genes:
            return genes, "three_algorithm_intersection"
    vf = BASE / "hub_genes_results" / "three_algorithm_vote.csv"
    vote = pd.read_csv(vf)
    genes = vote.loc[vote["vote"] >= 2, "gene_symbol"].tolist()[:max_n]
    return genes, "vote>=2_fallback"


def build_features(expr, sample, genes):
    cols = [g for g in genes if g in expr.index]
    sub = expr.loc[cols].T
    lab = (sample.loc[sub.index, "group"] == "PMOP").astype(int)
    return sub, lab.values, cols


def plot_roc_curves(roc_data, fname, title):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    colors = ["#A32D2D", "#185FA5", "#0F6E56", "#BA7517"]
    for i, (name, (fpr, tpr, roc_auc)) in enumerate(roc_data.items()):
        ax.plot(fpr, tpr, color=colors[i % len(colors)], lw=1.8,
                label=f"{name} (AUC={roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)


def nomogram(coef_df, fname):
    """文本化列线图（系数→score→概率），并输出可读的 SVG 风格条形图。"""
    fig, ax = plt.subplots(figsize=(8, 0.9 + 0.5 * len(coef_df)))
    genes = coef_df["gene_symbol"].tolist()
    coefs = coef_df["coefficient"].values
    ypos = np.arange(len(genes))
    ax.barh(ypos, np.abs(coefs), color=["#A32D2D" if c > 0 else "#185FA5" for c in coefs])
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes)
    ax.set_xlabel("|logistic coefficient| (red: risk, blue: protective)")
    ax.set_title("Nomogram - point contribution by gene")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)


def calibration_curve(y, prob, fname, n_bins=5):
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.digitize(prob, bins) - 1
    obs, pred = [], []
    for b in range(n_bins):
        m = idx == b
        if m.sum() > 0:
            pred.append(prob[m].mean())
            obs.append(y[m].mean())
        else:
            pred.append(np.nan); obs.append(np.nan)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.plot(pred, obs, "o-", color="#185FA5", lw=1.8)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed fraction")
    ax.set_title(f"Calibration (Brier={brier_score_loss(y, prob):.3f})")
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)
    return pred, obs


def dca_curve(y, prob, fname):
    thresholds = np.linspace(0.05, 0.95, 46)
    n = len(y)
    net_model, net_all, net_none = [], [], []
    for t in thresholds:
        pred = prob >= t
        tp = np.sum(pred & (y == 1)); fp = np.sum(pred & (y == 0))
        nb_model = tp / n - (fp / n) * (t / (1 - t))
        net_model.append(nb_model)
        net_all.append(np.mean(y) - (1 - np.mean(y)) * (t / (1 - t)))
        net_none.append(0.0)
    fig, ax = plt.subplots(figsize=(6, 4.8))
    ax.plot(thresholds, net_model, color="#A32D2D", lw=1.8, label="Model")
    ax.plot(thresholds, net_all, color="#185FA5", lw=1.4, ls="--", label="Treat all")
    ax.plot(thresholds, net_none, color="#888780", lw=1.2, label="Treat none")
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.set_ylim(-0.1, 0.4)
    ax.set_title("Decision curve analysis")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)


def external_roc(expr_ext, sample_ext, genes, coef_df, name):
    cols = [g for g in genes if g in expr_ext.index and g in set(coef_df["gene_symbol"])]
    if not cols:
        return None, 0, []
    X = expr_ext.loc[cols].T
    # 标准化参数必须与训练集一致——此处简化：用外部集自身 z-score（跨平台时可比性有限，已标注）
    X = (X - X.mean()) / X.std(ddof=0).replace(0, 1)
    cmap = coef_df.set_index("gene_symbol")["coefficient"].to_dict()
    score = sum(cmap[g] * X[g] for g in cols)
    y = (sample_ext.loc[X.index, "group"] == "PMOP").astype(int).values
    if len(np.unique(y)) < 2:
        return None, len(y), list(cols)
    prob = 1 / (1 + np.exp(-score))
    fpr, tpr, _ = roc_curve(y, prob)
    return (fpr, tpr, auc(fpr, tpr)), len(y), cols


def main():
    genes, src = load_hub_genes()
    print(f"Hub Genes ({src}): {genes}")
    expr = pd.read_csv(BASE / "GSE56815_exprs_matrix.csv", index_col=0)
    sample = pd.read_csv(BASE / "GSE56815_sample_info.csv", index_col=0)
    X, y, cols = build_features(expr, sample, genes)
    print(f"训练集特征: {X.shape[1]} genes, 样本 {len(y)}")

    pipe = Pipeline([("scaler", StandardScaler()),
                     ("lr", LogisticRegression(C=1.0, max_iter=5000, random_state=42))])
    cv_prob = cross_val_predict(pipe, X.values, y, cv=CV, method="predict_proba")[:, 1]
    fpr, tpr, _ = roc_curve(y, cv_prob)
    roc_auc = auc(fpr, tpr)
    roc_data = {"Training (5-fold CV)": (fpr, tpr, roc_auc)}

    # 全数据拟合模型系数用于列线图
    scaler = StandardScaler().fit(X.values)
    Xs = scaler.transform(X.values)
    lr = LogisticRegression(C=1.0, max_iter=5000, random_state=42).fit(Xs, y)
    coef_df = pd.DataFrame({"gene_symbol": cols, "coefficient": lr.coef_[0]})
    coef_df["odds_ratio"] = np.exp(coef_df["coefficient"])
    coef_df.to_csv(OUT / "model_coefficients.csv", index=False)
    print("\n模型系数:")
    print(coef_df.to_string(index=False))

    # 外部验证
    ext_summary = []
    g7429_e = BASE / "GSE7429_exprs_matrix.csv"
    if g7429_e.exists() and (BASE / "GSE7429_sample_info.csv").exists():
        ee = pd.read_csv(g7429_e, index_col=0)
        es = pd.read_csv(BASE / "GSE7429_sample_info.csv", index_col=0)
        roc, n, used = external_roc(ee, es, genes, coef_df, "GSE7429")
        if roc:
            roc_data["GSE7429 (external)"] = roc
            ext_summary.append({"dataset": "GSE7429", "n": n, "auc": roc[2], "genes_used": len(used)})
    g56814_e = BASE / "GSE56814_exprs_matrix.csv"
    if g56814_e.exists() and (BASE / "GSE56814_sample_info.csv").exists():
        ee = pd.read_csv(g56814_e, index_col=0)
        es = pd.read_csv(BASE / "GSE56814_sample_info.csv", index_col=0)
        roc, n, used = external_roc(ee, es, genes, coef_df, "GSE56814")
        if roc:
            roc_data["GSE56814 (cross-platform)"] = roc
            ext_summary.append({"dataset": "GSE56814", "n": n, "auc": roc[2], "genes_used": len(used)})

    plot_roc_curves(roc_data, "ROC_training_external.png",
                    "Diagnostic model ROC - training and external validation")
    nomogram(coef_df, "Nomogram.png")
    calibration_curve(y, cv_prob, "Calibration.png")
    dca_curve(y, cv_prob, "DCA.png")

    summary = {
        "hub_genes": genes, "hub_source": src,
        "training_cv_auc": roc_auc,
        "external": ext_summary,
        "note": "外部集标准化使用各自 z-score；跨平台 GSE56814 与同队列 GSE56815 非独立队列，仅作一致性参考。"
    }
    with open(OUT / "analysis_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    pd.DataFrame(ext_summary).to_csv(OUT / "external_validation.csv", index=False)
    print("\n汇总:")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
