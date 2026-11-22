#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第三阶段：LASSO 回归筛选 PMOP 核心特征基因（GSE56815，20 PMOP vs 20 Control）

方法：
  1) 特征集：
     - FDR_DEG333 : FDR<0.05 的 333 个 DEGs（主分析）
     - Nominal_P2581 : P<0.05 的 2581 个 DEGs（补充）
     - DRG_intersect : DRGs-DEGs 交集及扩展 DRG 中 P<0.05 的基因（机制聚焦）
  2) 标准化：每个基因 z-score（LASSO 对尺度敏感，必须标准化）
  3) 模型：
     - 主模型：L1 惩罚逻辑回归 LogisticRegressionCV（saga, StratifiedKFold=10）
       直接对二分类结局建模，等价于 glmnet family=binomial
     - 对照：LassoCV（gaussian，线性概率），等价于 glmnet family=gaussian
  4) 稳定性选择：RepeatedStratifiedKFold 下 100 次重采样（固定 C），统计每基因被选中频率
  5) 评估：交叉验证 AUC、ROC 曲线、系数路径图

输出：
  lasso_results/ 目录下的特征基因表、系数表、稳定性表、图与摘要
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

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Inconsistent values.*")

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "py_lib"))
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression, LassoCV, lasso_path
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.preprocessing import StandardScaler

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "lasso_results"
OUT.mkdir(exist_ok=True)

EXPR = BASE / "GSE56815_exprs_matrix.csv"
SAMPLE = BASE / "GSE56815_sample_info.csv"
DEG = BASE / "GSE56815_DEG_limma.csv"
GMT = BASE / "disulfidptosis_geneset_extended.gmt"

SEED = 42
N_SPLITS = 10
N_REPEATS = 100
C_VALS = np.logspace(-2, 2, 30)   # 逻辑回归正则强度倒数 C（越小正则越强）
LAMBDA_N = 100


def load_feature_sets(deg, expr_index):
    """构建多个候选特征集。"""
    fdr = deg[deg["adj.P.Val"] < 0.05].index
    nominal = deg[deg["P.Value"] < 0.05].index
    # DRG 聚焦集：扩展 53 DRG 中在 P<0.05 的基因（机制相关）
    drg = set()
    for line in open(GMT, encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "DISULFIDPTOSIS_EXTENDED":
            drg = set(parts[2:])
    drg_sig = sorted((drg & set(nominal)) & set(expr_index))
    return {
        "FDR_DEG333": list(fdr),
        "Nominal_P2581": list(nominal),
        "DRG_focus": drg_sig,
    }


def stability_selection(X, y, feature_names, best_C, n_repeats=N_REPEATS, seed=SEED):
    """稳定性选择：重复分层 CV，用固定 C 拟合，统计每个基因被选中（系数≠0）的频率。

    这是 stability selection 的标准做法：正则强度固定，仅重复数据划分，
    考察特征选择对样本扰动的稳健性，避免每次嵌套调参带来的偏差与开销。
    """
    counts = np.zeros(len(feature_names))
    cv = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats // 5,
                                      random_state=seed).split(X, y))
    for train_idx, _ in cv:
        model = LogisticRegression(C=best_C, penalty="l1", solver="saga",
                                   max_iter=5000, random_state=seed)
        model.fit(X[train_idx], y[train_idx])
        counts += (np.abs(model.coef_[0]) > 1e-8).astype(float)
    freq = counts / len(cv)
    return pd.DataFrame({"gene_symbol": feature_names, "selection_frequency": freq})


def run_logistic_l1(X, y, features, name):
    """L1 逻辑回归主模型 + 稳定性选择。"""
    print(f"\n=== [{name}] L1 逻辑回归 (n_features={len(features)}) ===")
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    model = LogisticRegressionCV(Cs=C_VALS, cv=cv, penalty="l1", solver="saga",
                                 scoring="roc_auc", max_iter=5000,
                                 random_state=SEED, n_jobs=1)
    model.fit(X, y)
    coef = model.coef_[0]
    chosen = np.abs(coef) > 1e-8
    best_C = float(model.C_[0])

    # 交叉验证预测 AUC（用选定 C 做 CV，避免嵌套偏差用简单 CV 近似）
    clf = LogisticRegression(C=best_C, penalty="l1", solver="saga",
                             max_iter=5000, random_state=SEED)
    y_score = cross_val_predict(clf, X, y, cv=cv, method="predict_proba", n_jobs=1)[:, 1]
    auc_cv = roc_auc_score(y, y_score)
    auc_train = roc_auc_score(y, model.predict_proba(X)[:, 1])

    coef_df = pd.DataFrame({
        "gene_symbol": features,
        "coefficient": coef,
        "abs_coefficient": np.abs(coef),
        "selected": chosen,
    }).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(OUT / f"LASSO_L1logistic_{name}_coefficients.csv", index=False)
    sel_df = coef_df[coef_df["selected"]].copy()

    stab = stability_selection(X, y, features, best_C)
    stab.to_csv(OUT / f"LASSO_{name}_stability.csv", index=False)
    merged = sel_df.merge(stab, on="gene_symbol", how="left")
    merged.to_csv(OUT / f"LASSO_{name}_selected_genes.csv", index=False)

    print(f"  best C = {best_C:.4f}")
    print(f"  选中基因数 = {int(chosen.sum())}")
    print(f"  CV AUC = {auc_cv:.4f} | train AUC = {auc_train:.4f}")
    if len(sel_df):
        print("  Top 系数基因:")
        print(merged.head(15)[["gene_symbol", "coefficient", "selection_frequency"]]
              .to_string(index=False))

    return {
        "name": name, "method": "L1-logistic", "n_features": len(features),
        "best_C": best_C, "n_selected": int(chosen.sum()),
        "cv_auc": auc_cv, "train_auc": auc_train,
        "selected_genes": ";".join(sel_df["gene_symbol"].tolist()),
    }, merged, y, y_score


def run_lasso_gaussian(X, y, features, name):
    """LassoCV（gaussian）经典对照。"""
    print(f"\n=== [{name}] LassoCV (gaussian) 对照 ===")
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    model = LassoCV(cv=cv, max_iter=50000,
                    random_state=SEED, n_jobs=1)
    model.fit(X, y)
    coef = model.coef_
    chosen = np.abs(coef) > 1e-8
    best_alpha = float(model.alpha_)

    coef_df = pd.DataFrame({
        "gene_symbol": features, "coefficient": coef,
        "abs_coefficient": np.abs(coef), "selected": chosen,
    }).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(OUT / f"LASSO_gaussian_{name}_coefficients.csv", index=False)
    sel = coef_df[coef_df["selected"]]

    print(f"  best alpha = {best_alpha:.4f}")
    print(f"  选中基因数 = {int(chosen.sum())}")
    if len(sel):
        print(sel.head(15)[["gene_symbol", "coefficient"]].to_string(index=False))
    return {
        "name": name, "method": "LassoCV-gaussian", "n_features": len(features),
        "best_alpha": best_alpha, "best_C": np.nan,
        "n_selected": int(chosen.sum()), "cv_auc": np.nan, "train_auc": np.nan,
        "selected_genes": ";".join(sel["gene_symbol"].tolist()),
    }


def plot_lasso_path(X, y, features, name, n_top=25):
    """绘制 LASSO 系数收缩路径（L1 逻辑回归，用 lasso_path 近似展示 gaussian 路径）。"""
    try:
        alphas, coefs, _ = lasso_path(X, y, n_alphas=60)
        # 取平均绝对系数最大的 top 基因画路径
        mean_abs = np.mean(np.abs(coefs), axis=1)
        idx = np.argsort(mean_abs)[::-1][:n_top]
        fig, ax = plt.subplots(figsize=(9, 5.5), dpi=150)
        for i in idx:
            ax.plot(alphas, coefs[i], lw=1.2, label=features[i])
        ax.set_xscale("log")
        ax.set_xlabel("alpha (regularization strength)")
        ax.set_ylabel("coefficient")
        ax.set_title(f"LASSO coefficient path — {name}")
        ax.legend(loc="upper right", fontsize=6, ncol=2)
        ax.grid(alpha=0.15)
        fig.tight_layout()
        fig.savefig(OUT / f"lasso_path_{name}.png")
        plt.close(fig)
    except Exception as e:
        print(f"[WARN] lasso path skipped: {e}")


def plot_roc(y, y_score, auc_cv, name):
    fpr, tpr, _ = roc_curve(y, y_score)
    fig, ax = plt.subplots(figsize=(5.5, 5), dpi=150)
    ax.plot(fpr, tpr, color="#D62728", lw=2, label=f"AUC = {auc_cv:.3f}")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"Cross-validated ROC — {name}")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(OUT / f"ROC_{name}.png")
    plt.close(fig)


def plot_stability(stab, name, top=25):
    s = stab.sort_values("selection_frequency", ascending=False).head(top)
    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    ax.barh(s["gene_symbol"][::-1], s["selection_frequency"][::-1], color="#1F77B4")
    ax.axvline(0.8, color="#D62728", ls="--", lw=1, label="80% selection")
    ax.set_xlabel("Selection frequency across 100 CV resamples")
    ax.set_title(f"LASSO stability selection — {name}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15, axis="x")
    fig.tight_layout()
    fig.savefig(OUT / f"stability_{name}.png")
    plt.close(fig)


def main():
    print("读取 GSE56815 数据 ...")
    expr = pd.read_csv(EXPR, index_col=0)
    info = pd.read_csv(SAMPLE, index_col=0).loc[list(expr.columns)]
    deg = pd.read_csv(DEG, index_col=0)
    y = (info["group"].to_numpy() == "PMOP").astype(int)
    print(f"  样本 {len(y)} 例 (PMOP {y.sum()} / Control {len(y)-y.sum()})")

    feature_sets = load_feature_sets(deg, expr.index)
    for k, v in feature_sets.items():
        print(f"  特征集 {k}: {len(v)} 基因")

    scaler = StandardScaler()
    X_all = expr.to_numpy(dtype=float).T   # samples x genes
    X_all = scaler.fit_transform(X_all)
    colmap = {g: i for i, g in enumerate(expr.index)}

    results = []
    # 主分析：FDR_DEG333
    primary_feats = feature_sets["FDR_DEG333"]
    idx = [colmap[g] for g in primary_feats if g in colmap]
    primary_feats = [g for g in primary_feats if g in colmap]
    Xp = X_all[:, idx]

    r, merged, yy, yscore = run_logistic_l1(Xp, y, primary_feats, "FDR_DEG333")
    results.append(r)
    plot_roc(yy, yscore, r["cv_auc"], "FDR_DEG333")
    stab_all = pd.read_csv(OUT / "LASSO_FDR_DEG333_stability.csv")
    plot_stability(stab_all, "FDR_DEG333")
    plot_lasso_path(Xp, y, primary_feats, "FDR_DEG333")
    results.append(run_lasso_gaussian(Xp, y, primary_feats, "FDR_DEG333"))

    # 补充：DRG_focus（机制聚焦）
    drg_feats = feature_sets["DRG_focus"]
    if len(drg_feats) >= 3:
        idx2 = [colmap[g] for g in drg_feats if g in colmap]
        drg_feats = [g for g in drg_feats if g in colmap]
        Xd = X_all[:, idx2]
        r_drg, _, _, _ = run_logistic_l1(Xd, y, drg_feats, "DRG_focus")
        results.append(r_drg)
        results.append(run_lasso_gaussian(Xd, y, drg_feats, "DRG_focus"))

    # 补充：Nominal_P2581（大特征集，仅 L1 逻辑回归）
    nom_feats = feature_sets["Nominal_P2581"]
    idx3 = [colmap[g] for g in nom_feats if g in colmap]
    nom_feats = [g for g in nom_feats if g in colmap]
    Xn = X_all[:, idx3]
    r_nom, _, _, _ = run_logistic_l1(Xn, y, nom_feats, "Nominal_P2581")
    results.append(r_nom)

    # 兜底：results 中若混入 (result_dict, ...) 形式的元组，取其中的 dict。
    # 否则 pd.DataFrame 内部会对元组执行 dict(tuple)，把首个 8 键的 dict
    # 当成键值对序列，抛出 "dictionary update sequence element #0 has length 8; 2 is required"。
    clean_results = []
    for rec in results:
        if isinstance(rec, dict):
            clean_results.append(rec)
        elif isinstance(rec, (tuple, list)) and len(rec) and isinstance(rec[0], dict):
            clean_results.append(rec[0])
        else:
            raise TypeError(f"无法识别的结果项类型: {type(rec)!r}")
    res_df = pd.DataFrame(clean_results)
    res_df.to_csv(OUT / "LASSO_summary.csv", index=False)

    # 核心基因：主分析中被选中且稳定性>=0.8 的基因
    core = merged[(merged["selected"]) & (merged["selection_frequency"] >= 0.8)]
    core = core.sort_values(["selection_frequency", "abs_coefficient"], ascending=False)
    core.to_csv(OUT / "core_feature_genes.csv", index=False)

    lines = [
        "第三阶段：LASSO 回归筛选 PMOP 核心特征基因（GSE56815）",
        "=" * 70,
        f"样本: {len(y)} 例 (PMOP {y.sum()} / Control {len(y)-y.sum()})，绝经后亚组",
        "特征标准化: 每个基因 z-score",
        "",
        "[模型]",
        "  主模型: L1 惩罚逻辑回归 (saga, 10-fold Stratified CV, 等价 glmnet binomial)",
        "  对照: LassoCV gaussian (线性概率, 等价 glmnet gaussian)",
        "  稳定性选择: 100 次重复分层 CV，统计基因被选中频率",
        "",
        "[各特征集结果]",
        res_df[["name", "method", "n_features", "n_selected", "cv_auc", "train_auc"]]
            .to_string(index=False),
        "",
        "[核心特征基因] (FDR_DEG333 主分析中 被选中 且 稳定性>=80%)",
        f"共 {len(core)} 个:",
        core[["gene_symbol", "coefficient", "selection_frequency"]].head(30).to_string(index=False)
        if len(core) else "  无基因同时满足两条件，见 LASSO_FDR_DEG333_selected_genes.csv",
        "",
        "[说明]",
        "  - 特征数远大于样本数（333 vs 40），LASSO 通过 L1 正则化实现特征压缩与选择。",
        "  - 稳定性选择用于降低单次数据划分的偶然性：频率越高，基因重要性越稳健。",
        "  - CV AUC 为 10 折分层交叉验证的预测性能，衡量特征组合的分类能力。",
    ]
    (OUT / "analysis_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n结果已保存到:", OUT)
    print(f"核心特征基因数 (稳定性>=80%): {len(core)}")
    if len(core):
        print(", ".join(core["gene_symbol"].tolist()))


if __name__ == "__main__":
    main()
