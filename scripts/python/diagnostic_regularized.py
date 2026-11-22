#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""第四阶段补充：强正则精简模型（缓解 18 基因/40 样本过拟合）。
对比 C=1.0（原）与 C=0.1/0.01（强 L2）在训练集 5-fold CV 与外部集的表现。
"""
import pathlib
import warnings
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "diagnostic_results"
import sys
sys.path.insert(0, str(BASE / "py_lib"))
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, auc

CV = StratifiedKFold(5, shuffle=True, random_state=42)
hub = pd.read_csv(BASE / "hub_genes_results" / "HubGenes_three_algorithm_intersection.csv")["gene_symbol"].tolist()
expr = pd.read_csv(BASE / "GSE56815_exprs_matrix.csv", index_col=0)
sample = pd.read_csv(BASE / "GSE56815_sample_info.csv", index_col=0)
X = expr.loc[hub].T.values
y = (sample.loc[expr.columns, "group"] == "PMOP").astype(int).values

# 外部集
ext = {}
for name in ["GSE7429", "GSE56814"]:
    ef = BASE / f"{name}_exprs_matrix.csv"
    sf = BASE / f"{name}_sample_info.csv"
    if ef.exists() and sf.exists():
        ee = pd.read_csv(ef, index_col=0)
        es = pd.read_csv(sf, index_col=0)
        genes = [g for g in hub if g in ee.index]
        Xe = ee.loc[genes].T
        ye = (es.loc[Xe.index, "group"] == "PMOP").astype(int).values
        ext[name] = (Xe.values, ye, genes)

rows = []
roc_data = {}
for C in [1.0, 0.1, 0.01]:
    pipe = Pipeline([("s", StandardScaler()), ("lr", LogisticRegression(C=C, max_iter=5000, random_state=42))])
    cvp = cross_val_predict(pipe, X, y, cv=CV, method="predict_proba")[:, 1]
    cv_auc = roc_auc_score(y, cvp)
    fpr, tpr, _ = roc_curve(y, cvp)
    roc_data[f"C={C} (CV)"] = (fpr, tpr, auc(fpr, tpr))
    row = {"C": C, "train_cv_auc": cv_auc}
    # 外部：全训练拟合后应用（各自 z-score）
    pipe.fit(X, y)
    sc = StandardScaler().fit(X)
    for ename, (Xe, ye, genes) in ext.items():
        # 按训练集 Hub 全列顺序对齐：缺失基因补 0，再用训练集 scaler 参数
        aligned = np.zeros((Xe.shape[0], len(hub)))
        for j, g in enumerate(hub):
            if g in genes:
                aligned[:, j] = Xe[:, genes.index(g)]
        Xs = (aligned - sc.mean_) / sc.scale_
        prob = pipe.predict_proba(Xs)[:, 1]
        e_auc = roc_auc_score(ye, prob)
        row[f"{ename}_auc"] = e_auc
        fpr2, tpr2, _ = roc_curve(ye, prob)
        roc_data[f"C={C} {ename}"] = (fpr2, tpr2, auc(fpr2, tpr2))
    rows.append(row)

rdf = pd.DataFrame(rows)
rdf.to_csv(OUT / "regularization_comparison.csv", index=False)
print(rdf.to_string(index=False))

fig, ax = plt.subplots(figsize=(6.5, 5.5))
colors = ["#A32D2D", "#185FA5", "#0F6E56"]
for i, (name, (fpr, tpr, a)) in enumerate(roc_data.items()):
    ax.plot(fpr, tpr, color=colors[i % 3], lw=1.5, ls="-" if "CV" in name else "--",
            label=f"{name} (AUC={a:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1)
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_title("Regularization effect on overfitting")
ax.legend(fontsize=7, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "ROC_regularization_comparison.png", dpi=140)
plt.close(fig)
print("done")
