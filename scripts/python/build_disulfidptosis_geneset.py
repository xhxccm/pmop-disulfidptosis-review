"""
构建双硫死亡（disulfidptosis）相关基因集，并与 GSE56814 表达矩阵取交集。

来源：
1) Liu X. et al., Nat Cell Biol 2023, PMID 36747082 —— 28 个 DRGs（原始实验证据 + CRISPR 筛选）
   经 PMID 36747082 来源文献整理（见 Comprehensive analysis ... ulcerative colitis, PMC11231342）
2) Li J. et al., Sci Rep 2023, DOI 10.1038/s41598-023-40381-5 —— 29 个 DRGs（上述 28 个 + RAC1）
3) CORE10 —— 上述 28 个中的核心机制子集（胱氨酸转运 + WAVE 复合体 + 肌动蛋白骨架）

输出：
- disulfidptosis_geneset.gmt
- disulfidptosis_geneset.csv
- GSE56814_DRG_expression.csv
- GSE56814_DRG_summary.csv
- GSE56814_DRG_score.csv（探索性：可用 DRG 的 z-score 均值）
"""

import pathlib

import pandas as pd

BASE = pathlib.Path(__file__).resolve().parent

DRG_LIU2023_28 = [
    "SLC7A11", "SLC3A2", "NCKAP1", "WASF2", "CYFIP1", "ABI2", "BRK1",
    "NUBPL", "NDUFA11", "LRPPRC", "OXSM", "NDUFS1", "GYS1", "FLNA",
    "FLNB", "MYL6", "MYH9", "MYH10", "ACTB", "ACTN4", "CAPZB", "CD2AP",
    "DSTN", "TLN1", "INF2", "PDLIM1", "IQGAP1", "RPN1",
]

DRG_LI2023_29 = DRG_LIU2023_28 + ["RAC1"]

DRG_CORE10 = [
    "SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "CYFIP1",
    "ABI2", "BRK1", "WASF2", "ACTB", "ACTN4",
]

SETS = {
    "DISULFIDPTOSIS_LIU2023_28": (
        "Liu X et al. Nat Cell Biol 2023 PMID:36747082; "
        "28 disulfidptosis-related genes (curated from PMID:36747082 / PMC11231342)",
        DRG_LIU2023_28,
    ),
    "DISULFIDPTOSIS_LI2023_29": (
        "Li J et al. Sci Rep 2023 DOI:10.1038/s41598-023-40381-5; "
        "29 DRGs (Liu2023 28 + RAC1)",
        DRG_LI2023_29,
    ),
    "DISULFIDPTOSIS_CORE10": (
        "Mechanistic core subset of Liu2023 DRGs: cystine transport + WRC + actin cytoskeleton",
        DRG_CORE10,
    ),
}


def write_gmt(path: pathlib.Path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, (desc, genes) in SETS.items():
            fh.write("\t".join([name, desc] + genes) + "\n")


def write_csv(path: pathlib.Path) -> pd.DataFrame:
    rows = []
    for name, (_, genes) in SETS.items():
        for g in genes:
            rows.append({"gene_set": name, "gene_symbol": g})
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df


def main():
    print("[1/4] 写出基因集文件 ...")
    write_gmt(BASE / "disulfidptosis_geneset.gmt")
    geneset_df = write_csv(BASE / "disulfidptosis_geneset.csv")
    for name, (_, genes) in SETS.items():
        print(f"      {name}: {len(genes)} genes")

    print("[2/4] 读取 GSE56814 表达矩阵 ...")
    expr = pd.read_csv(BASE / "GSE56814_exprs_matrix.csv", index_col=0)
    available = set(expr.index.astype(str))

    print("[3/4] 交集与覆盖度汇总 ...")
    summary_rows = []
    for name, (_, genes) in SETS.items():
        present = [g for g in genes if g in available]
        missing = [g for g in genes if g not in available]
        summary_rows.append({
            "gene_set": name,
            "n_genes": len(genes),
            "n_present_in_GSE56814": len(present),
            "coverage": round(len(present) / len(genes), 4),
            "present_genes": ";".join(present),
            "missing_genes": ";".join(missing),
        })
        print(f"      {name}: {len(present)}/{len(genes)} 存在于表达矩阵")
        if missing:
            print(f"        缺失: {', '.join(missing)}")
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(BASE / "GSE56814_DRG_summary.csv", index=False)

    # 用主基因集 28 个提取交集表达矩阵
    primary_present = [g for g in DRG_LIU2023_28 if g in available]
    drg_expr = expr.loc[primary_present].copy()
    drg_expr.to_csv(BASE / "GSE56814_DRG_expression.csv")

    print("[4/4] 计算探索性 DRG 评分（z-score 均值） ...")
    sample_info = pd.read_csv(BASE / "GSE56814_sample_info.csv", index_col=0)
    z = drg_expr.sub(drg_expr.mean(axis=1), axis=0).div(drg_expr.std(axis=1), axis=0)
    score = pd.DataFrame({
        "sample": z.columns,
        "DRG_score_zmean": z.mean(axis=0).values,
    }).set_index("sample")
    score = score.join(sample_info[["group"]])
    score.to_csv(BASE / "GSE56814_DRG_score.csv")
    print(score.groupby("group")["DRG_score_zmean"].agg(["count", "mean", "std"]).round(4))

    print("完成。")


if __name__ == "__main__":
    main()
