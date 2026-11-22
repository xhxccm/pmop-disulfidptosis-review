"""从已产出的 LASSO 结果文件生成汇总表与报告（不重跑模型）。"""
import pathlib
import pandas as pd

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "lasso_results"

# 各模型关键参数与性能（来自 gse56815_lasso.py 实际运行日志）
model_stats = [
    {"dataset": "GSE56815", "feature_set": "FDR_DEG333", "n_features": 333,
     "method": "L1-logistic", "best_C": 14.8735, "best_alpha": None,
     "n_selected": 101, "cv_auc": 0.9675, "train_auc": 1.0000},
    {"dataset": "GSE56815", "feature_set": "FDR_DEG333", "n_features": 333,
     "method": "LassoCV-gaussian", "best_C": None, "best_alpha": None,
     "n_selected": None, "cv_auc": None, "train_auc": None},
    {"dataset": "GSE56815", "feature_set": "DRG_focus", "n_features": 9,
     "method": "L1-logistic", "best_C": 0.2395, "best_alpha": None,
     "n_selected": 3, "cv_auc": 0.8350, "train_auc": 0.9275},
    {"dataset": "GSE56815", "feature_set": "DRG_focus", "n_features": 9,
     "method": "LassoCV-gaussian", "best_C": None, "best_alpha": 0.0366,
     "n_selected": 6, "cv_auc": None, "train_auc": None},
    {"dataset": "GSE56815", "feature_set": "Nominal_P2581", "n_features": 2581,
     "method": "L1-logistic", "best_C": 0.1743, "best_alpha": None,
     "n_selected": 15, "cv_auc": 0.8650, "train_auc": 1.0000},
]
# gaussian FDR 的选中数从系数文件读
g_fdr = pd.read_csv(OUT / "LASSO_gaussian_FDR_DEG333_coefficients.csv")
model_stats[1]["best_alpha"] = None
model_stats[1]["n_selected"] = int(g_fdr["selected"].sum())

sumdf = pd.DataFrame(model_stats)
sumdf.to_csv(OUT / "LASSO_summary.csv", index=False)

# 核心特征基因：主分析中被选中且稳定性>=80%
sel = pd.read_csv(OUT / "LASSO_FDR_DEG333_selected_genes.csv")
core = sel[(sel["selected"]) & (sel["selection_frequency"] >= 0.8)].copy()
core = core.sort_values(["selection_frequency", "abs_coefficient"], ascending=False)
core.to_csv(OUT / "core_feature_genes.csv", index=False)

# DRG 焦点集核心基因
drg_core = pd.read_csv(OUT / "LASSO_DRG_focus_selected_genes.csv")
drg_core = drg_core[drg_core["selected"]].sort_values("abs_coefficient", ascending=False)

lines = [
    "第三阶段：LASSO 回归筛选 PMOP 核心特征基因（GSE56815）",
    "=" * 72,
    "样本: 40 例 (PMOP 20 / Control 20)，绝经后亚组",
    "特征标准化: 每个基因 z-score",
    "",
    "[方法]",
    "  主模型: L1 惩罚逻辑回归 (saga, 10-fold Stratified CV, 等价 glmnet family=binomial)",
    "  对照: LassoCV gaussian (线性概率, 等价 glmnet family=gaussian)",
    "  稳定性选择: 100 次重复分层 CV（固定 C），统计每基因被选中频率",
    "",
    "[各特征集 / 模型结果]",
    sumdf[["feature_set", "method", "n_features", "n_selected", "cv_auc", "train_auc"]]
        .to_string(index=False),
    "",
    "[核心特征基因] 主分析 FDR_DEG333 中 被选中 且 稳定性 >= 80%",
    f"共 {len(core)} 个:",
    core[["gene_symbol", "coefficient", "selection_frequency"]].to_string(index=False),
    "",
    "[机制聚焦] DRG_focus 特征集（9 个 DRG-DEGs 交集基因）L1 筛选结果",
    "  CV AUC = 0.8350；选中 3 个基因:",
    drg_core[["gene_symbol", "coefficient", "selection_frequency"]].to_string(index=False),
    "",
    "[解读]",
    "  - FDR_DEG333 主模型 CV AUC=0.9675，说明 333 个 DEGs 的组合对 PMOP 有很强的分类判别力；",
    "    但 train AUC=1.0 提示存在过拟合（特征数 333 远大于样本数 40），需外部验证。",
    "  - 稳定性选择给出稳健核心：BAZ2A、MEX3D、MAP3K1、FLT3、PRDM14、TBC1D2 等",
    "    选择频率 >=0.93，是重复采样下最稳定的特征。",
    "  - DRG 聚焦集仅用 9 个交集基因即达 CV AUC=0.835，其中 PGRMC2 与 FLNC",
    "    选择频率 0.95，与差异分析/交集筛选结果一致，是连接双硫死亡机制的核心基因。",
    "  - Nominal_P2581 大特征集 CV AUC=0.865，但选中基因数更少(15)且频率偏低，",
    "    提示特征空间过大时信号被稀释。",
    "",
    "[交付]",
    "  LASSO_summary.csv / core_feature_genes.csv / LASSO_*_coefficients.csv",
    "  LASSO_*_stability.csv / LASSO_*_selected_genes.csv",
    "  ROC_FDR_DEG333.png / lasso_path_FDR_DEG333.png / stability_FDR_DEG333.png",
    "  gse56815_lasso.py（可复现）/ GSE56815_LASSO_glmnet.R（本地 R 复核）",
]
(OUT / "analysis_summary.txt").write_text("\n".join(lines), encoding="utf-8")
print("汇总完成")
print("核心特征基因数:", len(core))
print(", ".join(core["gene_symbol"].tolist()))
print("\n模型结果表:")
print(sumdf[["feature_set", "method", "n_selected", "cv_auc"]].to_string(index=False))
