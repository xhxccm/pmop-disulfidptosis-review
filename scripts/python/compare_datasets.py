"""定量对比 GSE56814 与 GSE56815 的差异信号强度，并检查 DRGs-DEGs 交集。"""
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

BASE = pathlib.Path(__file__).resolve().parent

new = pd.read_csv(BASE / "GSE56815_DEG_limma.csv", index_col=0)
old = pd.read_csv(BASE / "GSE56814_DEG_limma.csv", index_col=0)


def load_gmt(path):
    d = {}
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        d[f[0]] = set(f[2:])
    return d


gmt = load_gmt(BASE / "disulfidptosis_geneset_extended.gmt")
print("=== 基因集大小 ===", {k: len(v) for k, v in gmt.items()})
print()

print("=== 新旧数据集信号对比 ===")
rows = []
for name, d in [("GSE56814(旧)", old), ("GSE56815(新)", new)]:
    p = d["P.Value"].to_numpy()
    fdr = d["adj.P.Val"].to_numpy()
    n05 = int((p < 0.05).sum())
    n01 = int((p < 0.01).sum())
    nfdr = int((fdr < 0.05).sum())
    rows.append({
        "dataset": name, "n_genes": len(d),
        "P<0.05": n05, "P<0.01": n01,
        "随机期望(P<0.05)": int(0.05 * len(p)),
        "FDR<0.05": nfdr,
        "min_adjP": f"{fdr.min():.3g}",
        "max|logFC|": f"{d['logFC'].abs().max():.3f}",
        "mean_P": f"{p.mean():.4f}",
    })
print(pd.DataFrame(rows).to_string(index=False))
print()

print("=== DRGs-DEGs 交集（GSE56815 各阈值）===")
thr_defs = {
    "FDR<0.05": new[new["adj.P.Val"] < 0.05],
    "P<0.05": new[new["P.Value"] < 0.05],
    "P<0.05&|logFC|>=0.3": new[(new["P.Value"] < 0.05) & (new["logFC"].abs() >= 0.3)],
}
for thr, d in thr_defs.items():
    print(f"[{thr}] n_DEG={len(d)}")
    for k, gs in gmt.items():
        inter = sorted(set(d.index) & gs)
        print(f"   {k}: {len(inter)} {inter}")
print()

print("=== 扩展 DRG 在新数据集的差异统计（按 P 排序 top15）===")
k = max(gmt, key=lambda x: len(gmt[x]))
gs = gmt[k]
print("使用集合:", k, "n=", len(gs))
d = new.loc[sorted(gs & set(new.index))].sort_values("P.Value")
print(d[["logFC", "t", "P.Value", "adj.P.Val"]].head(15).round(5).to_string())
print(f"DRG 中 P<0.05: {(d['P.Value'] < 0.05).sum()}/{len(d)}  "
      f"FDR<0.05: {(d['adj.P.Val'] < 0.05).sum()}")
print()

print("=== 样本间相关性诊断 ===")
for label, stem in [("旧 GSE56814", "GSE56814"), ("新 GSE56815", "GSE56815")]:
    e = pd.read_csv(BASE / f"{stem}_exprs_matrix.csv", index_col=0)
    C = np.corrcoef(e.to_numpy(dtype=float).T)
    np.fill_diagonal(C, np.nan)
    print(f"{label}: 样本间相关性 mean={np.nanmean(C):.4f} median={np.nanmedian(C):.4f} "
          f"({e.shape[1]} samples)")
print()

print("=== 标签置换检验（普通 t，20 vs 20 / 26 vs 16）===")
for label, stem in [("旧 GSE56814", "GSE56814"), ("新 GSE56815", "GSE56815")]:
    e = pd.read_csv(BASE / f"{stem}_exprs_matrix.csv", index_col=0)
    s = pd.read_csv(BASE / f"{stem}_sample_info.csv", index_col=0).loc[list(e.columns)]
    g = s["group"].to_numpy()
    X = e.to_numpy(dtype=float)
    rng = np.random.default_rng(0)
    cnt = []
    for _ in range(50):
        gp = rng.permutation(g)
        _, pv = stats.ttest_ind(X[:, gp == "PMOP"], X[:, gp == "Control"],
                                axis=1, equal_var=True)
        cnt.append((pv < 0.05).sum())
    _, pv_real = stats.ttest_ind(X[:, g == "PMOP"], X[:, g == "Control"],
                                 axis=1, equal_var=True)
    print(f"{label}: 真实标签 P<0.05={(pv_real < 0.05).sum()}, "
          f"置换均值={np.mean(cnt):.0f} [min {min(cnt)}, max {max(cnt)}]")
