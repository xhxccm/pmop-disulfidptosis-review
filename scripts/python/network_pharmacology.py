#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第六阶段 6.1 + 6.4：蒙药蓝刺头网络药理学

6.1 活性成分-双硫死亡靶点网络：成分 ->靶点 -> DRG / Hub Gene -> PMOP
6.4 成分-靶点-通路调控网络可视化
数据来源：mongolian_lancitou_compound_target.csv / compounds.csv
"""
from __future__ import annotations

import pathlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "network_results"
OUT.mkdir(exist_ok=True)
CT = BASE / "mongolian_lancitou_compound_target.csv"
COMP = BASE / "mongolian_lancitou_compounds.csv"
GMT = BASE / "disulfidptosis_geneset_extended.gmt"
DEG = BASE / "GSE56815_DEG_limma.csv"


def load_sets():
    drg = set()
    for line in open(GMT, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if p[0] == "DISULFIDPTOSIS_EXTENDED":
            drg = set(p[2:])
    deg = pd.read_csv(DEG, index_col=0)
    deg_fdr = set(deg.index[deg["adj.P.Val"] < 0.05])
    deg_nom = set(deg.index[deg["P.Value"] < 0.05])
    hubf = BASE / "hub_genes_results" / "HubGenes_three_algorithm_intersection.csv"
    if hubf.exists():
        hub = set(pd.read_csv(hubf)["gene_symbol"])
    else:
        hub = set(pd.read_csv(BASE / "hub_genes_results" / "three_algorithm_vote.csv").head(10)["gene_symbol"])
    return drg, deg_fdr, deg_nom, hub


def main():
    ct = pd.read_csv(CT)
    comp = pd.read_csv(COMP)
    drg, deg_fdr, deg_nom, hub = load_sets()
    print(f"DRG {len(drg)}, FDR-DEG {len(deg_fdr)}, Hub {len(hub)}")

    # 单体成分（有 SMILES 的）
    monomers = set(comp.loc[comp["smiles"].notna(), "name_en"])
    ct_mono = ct[ct["compound_en"].isin(monomers)].copy()
    print(f"单体成分: {sorted(monomers)}")
    print(f"单体成分-靶点关系: {len(ct_mono)}")

    # 靶点交集分层
    ct_all_targets = set(ct["target_symbol"])
    overlap_drg = ct_all_targets & drg
    overlap_hub = ct_all_targets & hub
    overlap_deg = ct_all_targets & deg_fdr
    print(f"成分靶点∩DRG: {sorted(overlap_drg)}")
    print(f"成分靶点∩Hub: {sorted(overlap_hub)}")
    print(f"成分靶点∩FDR-DEG: {sorted(overlap_deg)}")

    # 核心网络：单体成分 -> (DRG∪Hub∪FDR-DEG) 靶点
    core_targets = overlap_drg | overlap_hub | overlap_deg
    core = ct_mono[ct_mono["target_symbol"].isin(core_targets)]
    core_edges = core[["compound_en", "target_symbol", "evidence_level"]].drop_duplicates()
    core_edges.to_csv(OUT / "compound_target_core_edges.csv", index=False)

    # 成分-靶点-类别标注
    def annotate(t):
        tags = []
        if t in drg: tags.append("DRG")
        if t in hub: tags.append("Hub")
        if t in deg_fdr: tags.append("FDR-DEG")
        if t in deg_nom: tags.append("Nominal-DEG")
        return "+".join(tags) if tags else "other"

    core_edges["target_class"] = core_edges["target_symbol"].map(annotate)
    core_edges.to_csv(OUT / "compound_target_core_edges.csv", index=False)

    # 网络图（二分图：成分-靶点）
    plot_bipartite(core_edges, "Compound_target_network.png")

    # 靶点类别统计
    class_stat = core_edges["target_class"].value_counts()
    class_stat.to_csv(OUT / "target_class_statistics.csv", header=["count"])

    # 成分-靶点-通路（用 DEGs 富集到的 KEGG 通路作为功能注释）
    kegg_map = build_kegg_annotation(core_targets)
    ctp = []
    for _, r in core_edges.iterrows():
        for pw in kegg_map.get(r["target_symbol"], []):
            ctp.append({"compound": r["compound_en"], "target": r["target_symbol"],
                        "pathway": pw, "target_class": r["target_class"]})
    ctp_df = pd.DataFrame(ctp)
    if not ctp_df.empty:
        ctp_df.to_csv(OUT / "compound_target_pathway_network.csv", index=False)
        plot_ctp(ctp_df, "Compound_target_pathway_network.png")

    summary = {
        "monomers": sorted(monomers),
        "n_monomer_target_edges": len(core_edges),
        "overlap_DRG": sorted(overlap_drg),
        "overlap_Hub": sorted(overlap_hub),
        "overlap_FDR_DEG": sorted(overlap_deg),
        "n_core_targets": len(core_targets),
        "target_class_statistics": class_stat.to_dict(),
    }
    with open(OUT / "analysis_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n汇总:")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def build_kegg_annotation(targets):
    """从本地 KEGG.gmt 提取核心靶点参与的通路（取前3条/基因，避免图过密）。"""
    result = {}
    gmt = BASE / "genesets" / "KEGG.gmt"
    if not gmt.exists():
        return result
    for line in open(gmt, encoding="utf-8", errors="replace"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        pw = p[0].split("%")[-1].strip()
        genes = set(p[2:])
        for t in targets & genes:
            result.setdefault(t, []).append(pw)
    return {t: v[:3] for t, v in result.items()}


def plot_bipartite(edges, fname):
    comps = sorted(edges["compound_en"].unique())
    targets = sorted(edges["target_symbol"].unique())
    if not comps or not targets:
        return
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(targets))))
    cy = np.linspace(0.5, len(comps) - 0.5, len(comps))
    ty = np.linspace(0.5, len(targets) - 0.5, len(targets))
    for _, r in edges.iterrows():
        ci = comps.index(r["compound_en"])
        ti = targets.index(r["target_symbol"])
        color = {"DRG": "#A32D2D", "Hub": "#BA7517", "FDR-DEG": "#185FA5"}.get(
            r["target_class"].split("+")[0], "#888780")
        ax.plot([1, 3], [cy[ci], ty[ti]], color=color, lw=0.8, alpha=0.6)
    ax.scatter([1] * len(comps), cy, s=180, color="#0F6E56", zorder=3)
    ax.scatter([3] * len(targets), ty, s=90,
               color=["#A32D2D" if "DRG" in targets[i] else "#888780" for i in range(len(targets))], zorder=3)
    for i, c in enumerate(comps):
        ax.text(0.92, cy[i], c, ha="right", va="center", fontsize=9)
    for i, t in enumerate(targets):
        ax.text(3.08, ty[i], t, ha="left", va="center", fontsize=9)
    ax.set_xlim(0.3, 4.2); ax.set_ylim(-1, max(len(comps), len(targets)) + 0.5)
    ax.set_xticks([1, 3]); ax.set_xticklabels(["Compound", "Target"])
    ax.set_yticks([]); ax.set_title("Echinops compound - target network")
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)


def plot_ctp(df, fname):
    # 通路-靶点数量条形图
    stat = df.groupby("pathway").size().sort_values(ascending=False).head(15)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(stat.index[::-1], stat.values[::-1], color="#534AB7")
    ax.set_xlabel("Number of compound-target links")
    ax.set_title("Top KEGG pathways in compound-target-pathway network")
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
