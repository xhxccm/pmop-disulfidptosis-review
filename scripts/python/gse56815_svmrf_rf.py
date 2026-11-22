#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第三阶段 3.2-3.4：SVM-RFE + Random Forest + 三算法交集 Hub Genes（GSE56815）

输入特征集与 LASSO 保持一致：FDR_DEG333（主）/ DRG_focus（机制聚焦）
方法：
  3.2 SVM-RFE：LinearSVC + 递归特征消除，按 5-fold CV AUC 选最优特征子集
  3.3 Random Forest：Mean Decrease Gini + 置换重要性排序
  3.4 三算法交集：LASSO(stability>=0.8) ∩ SVM-RFE(top) ∩ RF(top, 双指标一致)
输出：hub_genes_results/ 目录
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
from sklearn.svm import LinearSVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "hub_genes_results"
OUT.mkdir(exist_ok=True)
EXPR = BASE / "GSE56815_exprs_matrix.csv"
SAMPLE = BASE / "GSE56815_sample_info.csv"
DEG = BASE / "GSE56815_DEG_limma.csv"
GMT = BASE / "disulfidptosis_geneset_extended.gmt"
LASSO_CORE = BASE / "lasso_results" / "core_feature_genes.csv"
SEED = 42
CV = StratifiedKFold(5, shuffle=True, random_state=SEED)


def load_data():
    expr = pd.read_csv(EXPR, index_col=0)
    sample = pd.read_csv(SAMPLE, index_col=0)
    deg = pd.read_csv(DEG, index_col=0)
    y = (sample.loc[expr.columns, "group"] == "PMOP").astype(int).values
    fdr = list(deg.index[deg["adj.P.Val"] < 0.05])
    fdr = [g for g in fdr if g in expr.index]
    drg = set()
    for line in open(GMT, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if p[0] == "DISULFIDPTOSIS_EXTENDED":
            drg = set(p[2:])
    nominal = set(deg.index[deg["P.Value"] < 0.05])
    drg_focus = sorted((drg & nominal) & set(expr.index))
    return expr, y, fdr, drg_focus


def scale(expr, genes):
    X = expr.loc[genes].T.values
    X = StandardScaler().fit_transform(X)
    return X


def svm_rfe(X, y, genes, min_feat=5, max_feat=120, step=2):
    """SVM-RFE：从全特征递归消除，每步评估 CV AUC，返回最优特征数及排名。"""
    n = X.shape[1]
    step = max(1, min(step, n // 10))
    candidates = sorted(set([min_feat] + list(range(min_feat, min(max_feat, n) + 1, step))))
    results = []
    ranking = None
    for k in candidates:
        svc = LinearSVC(C=1.0, random_state=SEED, max_iter=20000, dual="auto")
        rfe = RFE(svc, n_features_to_select=k, step=step)
        rfe.fit(X, y)
        auc_cv = cross_val_score(rfe, X, y, cv=CV, scoring="roc_auc", n_jobs=1).mean()
        results.append({"n_features": k, "cv_auc": auc_cv})
        ranking = rfe.ranking_
    res_df = pd.DataFrame(results)
    best = res_df.loc[res_df["cv_auc"].idxmax()]
    k_best = int(best["n_features"])
    # 用最优 k 重新拟合，得到最终选中基因与排名
    svc = LinearSVC(C=1.0, random_state=SEED, max_iter=20000, dual="auto")
    rfe = RFE(svc, n_features_to_select=k_best, step=step).fit(X, y)
    support = rfe.support_
    selected = [g for g, s in zip(genes, support) if s]
    rank_series = pd.Series(rfe.ranking_, index=genes).sort_values()
    res_df.to_csv(OUT / "SVMRFE_cv_auc_curve.csv", index=False)
    rank_series.to_csv(OUT / "SVMRFE_feature_ranking.csv", header=["ranking"])
    pd.DataFrame({"gene_symbol": selected}).to_csv(OUT / "SVMRFE_selected_genes.csv", index=False)
    return selected, rank_series, float(best["cv_auc"]), res_df


def random_forest(X, y, genes, top_n=60, seed=SEED):
    rf = RandomForestClassifier(n_estimators=1000, random_state=seed, n_jobs=1)
    rf.fit(X, y)
    gini = pd.Series(rf.feature_importances_, index=genes, name="gini_importance")
    # 置换重要性（训练后，5 次重复）
    perm = permutation_importance(rf, X, y, n_repeats=5, random_state=seed,
                                  scoring="roc_auc", n_jobs=1)
    permi = pd.Series(perm.importances_mean, index=genes, name="perm_importance")
    imp = pd.concat([gini, permi], axis=1)
    imp["gini_rank"] = imp["gini_importance"].rank(ascending=False, method="min")
    imp["perm_rank"] = imp["perm_importance"].rank(ascending=False, method="min")
    imp["rank_sum"] = imp["gini_rank"] + imp["perm_rank"]
    imp = imp.sort_values("rank_sum")
    imp.to_csv(OUT / "RF_feature_importance.csv")
    top = list(imp.index[:top_n])
    return top, imp


def plot_svm_curve(res_df, title, fname):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(res_df["n_features"], res_df["cv_auc"], marker="o", ms=3, color="#185FA5")
    best = res_df.loc[res_df["cv_auc"].idxmax()]
    ax.axvline(best["n_features"], color="#A32D2D", ls="--", lw=1)
    ax.set_xlabel("Number of features")
    ax.set_ylabel("CV AUC (5-fold)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=130)
    plt.close(fig)


def plot_rf_top(imp, title, fname, n=25):
    top = imp.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(top.index, top["gini_importance"], color="#0F6E56", alpha=0.85)
    ax.set_xlabel("Mean decrease Gini")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=130)
    plt.close(fig)


def run_set(expr, y, genes, tag):
    print(f"\n=== 特征集 {tag} ({len(genes)} genes) ===")
    X = scale(expr, genes)
    sel_svm, rank_svm, auc_svm, curve = svm_rfe(X, y, genes)
    print(f"  SVM-RFE: 最优 {len(sel_svm)} 基因, CV AUC={auc_svm:.4f}")
    top_rf, imp = random_forest(X, y, genes)
    print(f"  RF: top {len(top_rf)} 基因")
    plot_svm_curve(curve, f"SVM-RFE CV AUC ({tag})", f"SVMRFE_cv_auc_{tag}.png")
    plot_rf_top(imp, f"RF importance ({tag})", f"RF_importance_{tag}.png")
    return {"svm": sel_svm, "rf": top_rf, "svm_auc": auc_svm,
            "rank_svm": rank_svm, "imp": imp}


def main():
    expr, y, fdr, drg_focus = load_data()
    lasso = pd.read_csv(LASSO_CORE)
    lasso_genes = list(lasso.loc[lasso["selection_frequency"] >= 0.8, "gene_symbol"])
    print(f"LASSO 稳定核心基因 (>=0.8): {len(lasso_genes)}")

    results = {}
    for tag, genes in [("FDR_DEG333", fdr), ("DRG_focus", drg_focus)]:
        results[tag] = run_set(expr, y, genes, tag)

    # ---- 3.4 三算法交集（主分析 FDR_DEG333）----
    r = results["FDR_DEG333"]
    set_lasso = set(lasso_genes)
    set_svm = set(r["svm"])
    set_rf = set(r["rf"])
    inter3 = sorted(set_lasso & set_svm & set_rf)
    inter_ls = sorted(set_lasso & set_svm)
    inter_lr = sorted(set_lasso & set_rf)
    inter_sr = sorted(set_svm & set_rf)
    union = sorted(set_lasso | set_svm | set_rf)

    # 投票计分：出现在 n 个算法中即计 n 分
    vote = {}
    for g in union:
        vote[g] = (g in set_lasso) + (g in set_svm) + (g in set_rf)
    vote_df = pd.DataFrame({
        "gene_symbol": union,
        "vote": [vote[g] for g in union],
        "in_LASSO": [g in set_lasso for g in union],
        "in_SVMRFE": [g in set_svm for g in union],
        "in_RF": [g in set_rf for g in union],
    }).sort_values("vote", ascending=False)
    vote_df.to_csv(OUT / "three_algorithm_vote.csv", index=False)
    pd.DataFrame({"gene_symbol": inter3}).to_csv(OUT / "HubGenes_three_algorithm_intersection.csv", index=False)

    # 机制聚焦集交集
    rd = results["DRG_focus"]
    inter3_drg = sorted(set(lasso_genes) & set(rd["svm"]) & set(rd["rf"]))
    pd.DataFrame({"gene_symbol": inter3_drg}).to_csv(
        OUT / "HubGenes_DRGfocus_intersection.csv", index=False)

    summary = {
        "n_LASSO": len(set_lasso), "n_SVMRFE": len(set_svm), "n_RF": len(set_rf),
        "n_three_way": len(inter3), "three_way_genes": inter3,
        "n_LASSO_SVM": len(inter_ls), "LASSO_SVM": inter_ls,
        "n_LASSO_RF": len(inter_lr), "LASSO_RF": inter_lr,
        "n_SVM_RF": len(inter_sr), "SVM_RF": inter_sr,
        "DRGfocus_three_way": inter3_drg,
        "svm_best_cv_auc": r["svm_auc"],
    }
    with open(OUT / "analysis_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n=== 3.4 三算法交集 ===")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
