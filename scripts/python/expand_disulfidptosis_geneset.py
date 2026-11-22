#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
扩展双硫死亡（disulfidptosis）基因集至 50~100 个基因。

策略：汇总多个已发表来源的 DRG 列表，取并集，用 org.Hs.eg.db 规范符号，
剔除 NADPH 等非基因条目，按证据来源分层，构建扩展基因集。

来源（均为已发表文献，逐层累积去重）：
  S1  Liu X et al., Nat Cell Biol 2023, PMID 36747082 (CRISPR 筛选, 28 核心基因)
  S2  Li J et al., Sci Rep 2023, DOI 10.1038/s41598-023-40381-5 (28 + RAC1)
  S3  Zhang P et al., Aging 2024, PMID 38738994 (骨质疏松 RNA-seq, 7 诊断 DRG)
  S4  Pan-cancer analysis, Cancer Genet 2023, PMID 37879141 (15 骨架基因)
  S5  Front Cell Dev Biol 2026, DOI 10.3389/fcell.2026.1790068 (文献+GeneCards, 43)
  S6  PMC11596291 综述 Table 1 (23 基因, 含代谢/PPP 通路)
  S7  JIR 2024 骨关节炎研究 (GeneCards 33 DRG, 本文明确列出部分)
  S8  机制相关基因（低置信层）：Aging 2024 OP 与 2026 综述明确讨论的
      SLC7A11-TXNRD1/NFATc1 轴、PGRMC2、GLUT2/4、SLC26A2

输出基因集：
  DISULFIDPTOSIS_EXTENDED  — 全部来源并集（目标 50~100）
  DISULFIDPTOSIS_CORE10    — 核心机制子集（保留）
  DISULFIDPTOSIS_LIU2023_28 / LI2023_29 — 原集合（保留）
"""
from __future__ import annotations

import pathlib
import sqlite3

import pandas as pd

BASE = pathlib.Path(__file__).resolve().parent
ORG_SQL = BASE / "geo_cache" / "org.Hs.eg.sqlite"
DEG_CSV = BASE / "GSE56814_DEG_limma.csv"

# ---------- 各来源基因（文献原文） ----------
S1_LIU28 = """SLC7A11 SLC3A2 NCKAP1 WASF2 CYFIP1 ABI2 BRK1 NUBPL NDUFA11
LRPPRC OXSM NDUFS1 GYS1 FLNA FLNB MYL6 MYH9 MYH10 ACTB ACTN4 CAPZB
CD2AP DSTN TLN1 INF2 PDLIM1 IQGAP1 RPN1""".split()

S2_LI29 = S1_LIU28 + ["RAC1"]

S3_OP7 = "FLNA ACTB PRDX1 SLC7A11 NUBPL OXSM RAC1".split()

S4_PANC15 = "SLC7A11 INF2 CD2AP PDLIM1 ACTN4 MYH9 MYH10 IQGAP1 FLNA FLNB TLN1 MYL6 ACTB DSTN CAPZB".split()

S5_FRONTIERS43 = """CYFIP1 EMILIN1 AJAP1 TM9SF2 FLNA FLNB MYH9 TLN1 ACTB MYL6 MYH10
CAPZB DSTN IQGAP1 ACTN4 PDLIM1 CD2AP INF2 SLC7A11 GLUT NCKAP1 SLC3A2
GYS1 NDUFA11 RAC1 RPN1 NDUFS1 NUBPL LRPPRC OXSM NADPH BRK1 WASF2 ABI2
PRDX1 ACSL4 BAK1 TLN2 FLNC MYL6B DBN1 ACTN1 MYH11""".split()

S6_REVIEW23 = """SLC7A11 SLC3A2 GLUT1 GLUT3 NCKAP1 CYFIP1 WAVE2 ABI2 HSPC300 RAC1
RPN1 ATF4 NUBPL NDUFA11 LRPPRC OXSM NDUFS1 GYS1 G6PD 6-PGD TALDO1 TKT
PRC1""".split()

S7_OA = """SLC2A1 SLC2A3 OXSM GYS1 SLC3A2 RPN1 G6PD RAC1 CAPZB MYH9 MYL6
ACTB FLNA TM9SF2 DSTN SLC7A11""".split()

# 机制相关（Tier 3）：文献明确讨论的双硫死亡调控/效应基因，但非核心 CRISPR 筛选集
S8_MECH = """TXNRD1 NFATC1 PGRMC2 SLC26A2""".split()

SOURCES = {
    "S1_Liu2023_NCB": S1_LIU28,
    "S2_Li2023_SciRep": S2_LI29,
    "S3_Zhang2024_Aging_OP": S3_OP7,
    "S4_LiuTang2023_CancerGenet": S4_PANC15,
    "S5_Frontiers2026_MCD": S5_FRONTIERS43,
    "S6_Review_PMC11596291": S6_REVIEW23,
    "S7_OA2024_GeneCards": S7_OA,
    "S8_Mechanism_related": S8_MECH,
}

# 别名 -> 官方符号（基于 org.Hs.eg.db 验证；NADPH 是分子不是基因，剔除）
MANUAL_ALIAS = {
    "GLUT": "SLC2A1", "GLUT1": "SLC2A1", "GLUT2": "SLC2A2",
    "GLUT3": "SLC2A3", "GLUT4": "SLC2A4", "WAVE2": "WASF2",
    "HSPC300": "BRK1", "6-PGD": "PGD",
}
DROP = {"NADPH"}  # 还原型辅酶，非基因


def canonic(g: str, sym: set[str]) -> str | None:
    g = g.strip()
    if not g or g in DROP:
        return None
    if g in MANUAL_ALIAS:
        return MANUAL_ALIAS[g]
    return g if g in sym else None


def main() -> None:
    con = sqlite3.connect(ORG_SQL)
    sym = {r[0] for r in con.execute("SELECT symbol FROM gene_info WHERE symbol IS NOT NULL")}
    con.close()

    # 记录每个基因的来源与原始写法
    origin: dict[str, dict] = {}
    for src, genes in SOURCES.items():
        for g in genes:
            cg = canonic(g, sym)
            if cg is None:
                continue
            d = origin.setdefault(cg, {"gene": cg, "sources": [], "raw_names": set()})
            d["sources"].append(src)
            if g != cg:
                d["raw_names"].add(g)

    ext = sorted(origin)
    print(f"扩展基因集总数: {len(ext)}")
    assert 50 <= len(ext) <= 100, f"基因数 {len(ext)} 不在 50~100 区间，请核查来源"

    # 与表达矩阵覆盖
    deg = pd.read_csv(DEG_CSV, index_col=0)
    in_mat = [g for g in ext if g in deg.index]
    missing = [g for g in ext if g not in deg.index]
    print(f"在 GSE56814 表达矩阵中: {len(in_mat)}/{len(ext)}；缺失: {missing}")

    # 分层：核心（S1/S2 实验筛选）vs 扩展（仅文献/GeneCards）
    core = set(S2_LI29)
    mech = set(S8_MECH)
    def assign_tier(g: str) -> str:
        if g in core:
            return "T1_core_experimental"
        if g in mech:
            return "T3_mechanism_related"
        return "T2_literature_GeneCards"
    tier = {g: assign_tier(g) for g in ext}

    # --- 写出来源明细表 ---
    rows = []
    for g in ext:
        d = origin[g]
        rows.append({
            "gene_symbol": g,
            "tier": tier[g],
            "n_sources": len(set(d["sources"])),
            "sources": ";".join(sorted(set(d["sources"]))),
            "original_aliases": ";".join(sorted(d["raw_names"])),
            "in_GSE56814": g in deg.index,
        })
    detail = pd.DataFrame(rows).sort_values(["tier", "n_sources", "gene_symbol"],
                                           ascending=[True, False, True])
    detail.to_csv(BASE / "disulfidptosis_extended_evidence.csv", index=False)

    # --- 写出 GMT（第 2 字段为来源说明） ---
    def w_gmt(path: pathlib.Path, name: str, desc: str, genes: list[str], mode: str = "a") -> None:
        with open(path, mode, encoding="utf-8") as fh:
            fh.write(f"{name}\t{desc}\t" + "\t".join(genes) + "\n")

    (BASE / "disulfidptosis_geneset_extended.gmt").write_text("", encoding="utf-8")  # 清空
    w_gmt(BASE / "disulfidptosis_geneset_extended.gmt", "DISULFIDPTOSIS_EXTENDED",
          f"Extended DRG set: union of 7 published sources (Liu2023/Li2023/Aging2024/"
          f"CancerGenet2023/Frontiers2026/Review PMC11596291/OA2024), {len(ext)} genes, "
          f"symbols normalized by org.Hs.eg.db",
          ext)
    w_gmt(BASE / "disulfidptosis_geneset_extended.gmt", "DISULFIDPTOSIS_EXTENDED_IN_MATRIX",
          f"Extended DRG set restricted to genes detected in GSE56814 ({len(in_mat)} genes)",
          in_mat)
    # 追加保留原有三个集合
    for line in (BASE / "disulfidptosis_geneset.gmt").read_text(encoding="utf-8").splitlines():
        with open(BASE / "disulfidptosis_geneset_extended.gmt", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    # --- 写出 CSV 长表 ---
    csv_rows = []
    for name, genes in [("DISULFIDPTOSIS_EXTENDED", ext),
                        ("DISULFIDPTOSIS_EXTENDED_IN_MATRIX", in_mat)]:
        csv_rows += [{"gene_set": name, "gene_symbol": g} for g in genes]
    old = pd.read_csv(BASE / "disulfidptosis_geneset.csv")
    pd.concat([pd.DataFrame(csv_rows), old], ignore_index=True).to_csv(
        BASE / "disulfidptosis_geneset_extended.csv", index=False)

    # --- 摘要 ---
    src_counts = {s: sum(1 for g in ext if s in origin[g]["sources"]) for s in SOURCES}
    summary = [
        "扩展双硫死亡基因集构建摘要",
        "=" * 60,
        f"最终基因数: {len(ext)}（目标 50~100）",
        f"表达矩阵覆盖: {len(in_mat)}/{len(ext)}；缺失 {missing}",
        "",
        "各来源贡献（规范后基因数）:",
    ] + [f"  {k}: {v}" for k, v in src_counts.items()] + [
        "",
        f"证据分层: T1 实验筛选核心 {sum(1 for g in ext if tier[g].startswith('T1'))} 个；"
        f"T2 文献/GeneCards 扩展 {sum(1 for g in ext if tier[g].startswith('T2'))} 个；"
        f"T3 机制相关 {sum(1 for g in ext if tier[g].startswith('T3'))} 个",
        "",
        "符号规范化: GLUT/GLUT1→SLC2A1, GLUT3→SLC2A3, WAVE2→WASF2, "
        "HSPC300→BRK1, 6-PGD→PGD；NADPH（非基因）已剔除。",
        "",
        "完整基因列表:",
        "  " + " ".join(ext),
    ]
    (BASE / "disulfidptosis_extended_summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary[:16]))


if __name__ == "__main__":
    main()
