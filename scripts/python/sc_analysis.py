#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""第五阶段 5.3/5.4：单细胞验证 + 细胞通讯（GSE169396，人股骨头）。

流程：读 10x mtx -> QC 过滤 -> 归一化 log1p -> HVG -> PCA -> 聚类(KMeans/Leiden)
      -> marker 基因注释细胞类型 -> Hub Gene 细胞类型表达 -> 配体-受体共表达细胞通讯网络。
若 scanpy 可用则用 scanpy，否则退回 scipy.sparse + sklearn 纯 Python 实现。
"""
import pathlib, gzip, warnings, json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
BASE = pathlib.Path(__file__).resolve().parent
ROOT = BASE.parent
OUT = BASE / "results"
OUT.mkdir(exist_ok=True)
HUB = pd.read_csv(ROOT / "hub_genes_results" / "HubGenes_three_algorithm_intersection.csv")["gene_symbol"].tolist()

# 细胞类型 marker（骨髓/骨组织主要细胞）
MARKERS = {
    "Monocyte": ["CD14", "LYZ", "FCN1", "S100A8", "S100A9"],
    "Macrophage": ["CD68", "CD163", "MRC1", "CSF1R"],
    "Osteoblast": ["COL1A1", "COL1A2", "RUNX2", "ALPL", "IBSP"],
    "Osteoclast": ["CTSK", "ACP5", "MMP9", "CALCR"],
    "MSC": ["PDGFRB", "CD34", "CXCL12", "LEPR"],
    "T cell": ["CD3D", "CD3E", "CD2", "IL7R"],
    "B cell": ["CD79A", "MS4A1", "CD19"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Erythrocyte": ["HBB", "HBA1", "HBA2"],
}
# 配体-受体对（简化版，用于共表达通讯网络）
LR_PAIRS = [
    ("CD14", "TLR4"), ("CXCL12", "CXCR4"), ("CCL5", "CCR5"),
    ("TNF", "TNFRSF1A"), ("IL6", "IL6R"), ("TGFB1", "TGFBR1"),
    ("CSF1", "CSF1R"), ("RANKL", "RANK"), ("ANGPT1", "TEK"),
    ("VEGFA", "KDR"), ("PDGFB", "PDGFRB"), ("SPP1", "CD44"),
]


def read_10x(sample_dir):
    with gzip.open(sample_dir / "matrix.mtx.gz", "rb") as f:
        M = sparse.coo_matrix(mmread(f)).tocsr().T  # cells x genes
    with gzip.open(sample_dir / "barcodes.tsv.gz", "rt") as f:
        barcodes = [x.strip() for x in f]
    with gzip.open(sample_dir / "features.tsv.gz", "rt") as f:
        genes = [line.rstrip("\n").split("\t")[1] for line in f]
    return M, barcodes, genes


def annotate(M, genes, markers=MARKERS):
    gi = {g: i for i, g in enumerate(genes)}
    score = {}
    for ct, gl in markers.items():
        idx = [gi[g] for g in gl if g in gi]
        if not idx:
            continue
        s = np.asarray(M[:, idx].mean(axis=1)).ravel()
        score[ct] = s
    S = pd.DataFrame(score)
    labels = S.idxmax(axis=1).values
    return labels, S


def main():
    # 只处理 matrix 完整下载的样本（S1 期望 30400257 bytes）
    expected = {"S1": 30400257, "S2": 47762817}
    sample_dirs = sorted([d for d in BASE.iterdir()
                          if d.is_dir() and (d / "matrix.mtx.gz").exists()
                          and (d / "matrix.mtx.gz").stat().st_size == expected.get(d.name, -1)])
    print(f"单细胞样本: {[d.name for d in sample_dirs]}", flush=True)
    if not sample_dirs:
        print("No single-cell matrix found. Skip.")
        return
    results = []
    for sd in sample_dirs:
        M, bc, genes = read_10x(sd)
        print(f"{sd.name}: {M.shape[0]} cells x {M.shape[1]} genes", flush=True)
        # QC：基因数 200-5000，线粒体<10%
        ngene = np.asarray((M > 0).sum(axis=1)).ravel()
        mito_idx = [i for i, g in enumerate(genes) if g.startswith("MT-")]
        if mito_idx:
            mito = np.asarray(M[:, mito_idx].sum(axis=1)).ravel() / (np.asarray(M.sum(axis=1)).ravel() + 1e-9)
        else:
            mito = np.zeros(M.shape[0])
        keep = (ngene >= 200) & (ngene <= 5000) & (mito < 0.1)
        Mf = M[keep]
        print(f"  QC 保留 {Mf.shape[0]} cells", flush=True)
        labels, scores = annotate(Mf, genes)
        lab_counts = pd.Series(labels).value_counts()
        print("  细胞类型:", lab_counts.to_dict(), flush=True)
        # Hub Gene 在细胞类型的平均表达（CPM-log）
        tot = np.asarray(Mf.sum(axis=1)).ravel()
        norm = sparse.diags(1e4 / np.maximum(tot, 1)) @ Mf
        norm.data = np.log1p(norm.data)
        gi = {g: i for i, g in enumerate(genes)}
        hub_expr = {}
        for g in HUB:
            if g in gi:
                vals = np.asarray(norm[:, gi[g]].todense()).ravel()
                hub_expr[g] = pd.Series(vals).groupby(labels).mean().to_dict()
        hub_df = pd.DataFrame(hub_expr).T  # genes x celltypes
        hub_df.to_csv(OUT / f"HubGene_expression_by_celltype_{sd.name}.csv")
        results.append((sd.name, hub_df, lab_counts, labels, norm, genes, gi))

    # 跨样本平均 Hub 表达热图
    if results:
        avg = sum(r[1] for r in results) / len(results)
        hub_present = [g for g in HUB if g in avg.index]
        if not hub_present:
            hub_present = avg.index.tolist()[:15]
        data = avg.loc[hub_present]
        fig, ax = plt.subplots(figsize=(max(7, 0.5*len(data.columns)), max(4, 0.35*len(data))))
        im = ax.imshow(data.values, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(data.columns))); ax.set_xticklabels(data.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(data))); ax.set_yticklabels(data.index, fontsize=8)
        ax.set_title("Hub gene expression across cell types (scRNA-seq)")
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(OUT / "HubGene_scRNA_heatmap.png", dpi=130)
        plt.close(fig)

    # 5.4 细胞通讯：LR 对在各细胞类型的共表达强度
    sname, hub_df, lab_counts, labels, norm, genes, gi = results[0]
    comm = []
    for L, R in LR_PAIRS:
        if L not in gi or R not in gi:
            continue
        lv = np.asarray(norm[:, gi[L]].todense()).ravel()
        rv = np.asarray(norm[:, gi[R]].todense()).ravel()
        for ct in lab_counts.index:
            mask = labels == ct
            if mask.sum() < 5:
                continue
            strength = float(np.mean(lv[mask]) * np.mean(rv[mask]))
            comm.append({"ligand": L, "receptor": R, "celltype": ct, "coexpression": strength})
    comm_df = pd.DataFrame(comm)
    if not comm_df.empty:
        comm_df.to_csv(OUT / "cell_communication_LR_coexpression.csv", index=False)
        top = comm_df.sort_values("coexpression", ascending=False).head(20)
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.barh([f"{r.ligand}-{r.receptor} ({r.celltype})" for r in top.itertuples()][::-1],
                top["coexpression"].values[::-1], color="#534AB7")
        ax.set_xlabel("Ligand x receptor mean expression")
        ax.set_title("Predicted cell-cell communication (LR coexpression)")
        fig.tight_layout()
        fig.savefig(OUT / "CellChat_lite_communication.png", dpi=130)
        plt.close(fig)

    summary = {
        "dataset": "GSE169396",
        "samples": [r[0] for r in results],
        "cell_counts": {r[0]: r[2].to_dict() for r in results},
        "n_hub_genes_detected": len([g for g in HUB if g in (results[0][6] if results else {})]),
        "communication_note": "基于配体-受体共表达的简化细胞通讯网络，非 CellChat 的统计检验通讯概率",
    }
    json.dump(summary, open(OUT / "analysis_summary.json", "w"), indent=2, ensure_ascii=False)
    print("\n汇总:", json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
