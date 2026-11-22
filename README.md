# Reproducibility repository

> Companion code and data for the methodological review:
> "Pipelines, Pitfalls, and External Validation in Multi-omics Biomarker Discovery: A Worked Example from Disulfidptosis-Related Genes in Postmenopausal Osteoporosis"
>
> 本仓库收录综述案例研究的全部可复现脚本、数据与图表，供投稿时"Data and code availability"声明引用。

---

## 1. 仓库结构（建议上传布局）

```
repository/
├── README.md                      # 本文件
├── LICENSE                        # 建议 MIT 或 CC-BY-4.0
├── data/
│   ├── GSE56815_exprs_matrix.csv  # 主队列表达矩阵（40例，13,038基因）
│   ├── GSE56815_sample_info.csv   # 分组信息
│   ├── GSE7429_exprs_matrix.csv   # 外部验证1（B细胞）
│   ├── GSE7429_sample_info.csv
│   ├── GSE56814_exprs_matrix.csv  # 外部验证2（跨平台）
│   ├── GSE56814_sample_info.csv
│   ├── GSE56815_DEG_limma.csv     # limma 差异分析结果
│   ├── disulfidptosis_geneset_extended.gmt
│   └── diagnostic_results/        # 模型系数、外部验证、正则化对比
├── scripts/
│   ├── python/                    # 数据预处理、limma复现、机器学习、富集、网络药理
│   └── R/                         # R_native/reproduce_*.R（6个脚本）
├── figures/
│   ├── main/                      # 正文主图 Fig 1/2/3/4/5/6/7/8 (PNG/PDF/SVG)
│   └── supplementary/             # 90 项补充材料（见 supplementary_index.csv）
├── supplementary_index.csv        # 补充材料索引
└── environment/
    ├── requirements.txt           # Python 依赖及版本
    └── sessionInfo.txt            # R 会话信息（R 4.4.1 + 包版本）
```

> 注意：GEO 原始数据受数据库使用条款约束，本仓库仅存放**衍生的表达矩阵与结果表**。
> 上传前请确认各 GEO 数据集的许可允许再分发；若不允许，改为在 README 中给出 GEO 登录号与下载脚本。

---

## 2. 核心可复现结果

| 项目 | 结果 | 脚本 |
|---|---|---|
| limma 差异分析 | FDR<0.05: 333 基因 | `scripts/python/gse56814_de_limma.py`（对应 GSE56815）|
| DRG-DEG 交集 | 严格阈值 2 个 / 宽松阈值 9 个 | `build_drg_expression_GSE56815.py` |
| LASSO | 10折 CV AUC=0.9675；38 稳定性基因 | `gse56815_lasso.py` |
| SVM-RFE / RF | 67 / 60 基因 | `gse56815_svmrf_rf.py` |
| 三算法交集 | 18 Hub 基因 | `gse56815_svmrf_rf.py` |
| 诊断模型 | 训练 AUC=1.000；外部 0.340 / 0.572 | `gse56815_diagnostic_model.py` |
| 正则化对比 | C=1/0.1/0.01 外部 AUC 稳定 0.35–0.50 | `diagnostic_regularized.py` |
| R 交叉验证 | 三个 AUC 逐位一致 | `scripts/R/reproduce_fig4.R` |
| GSEA | fgsea 四库，KEGG 与 Python 近乎逐位吻合 | `R_native/reproduce_gsea_fgsea.R` |

## 3. 环境与版本（可复现性关键）

**Python**：3.13.12/3.13.14。精确版本清单见 **`full_requirements_exact.txt`（64 个包，从 `py_lib/*.dist-info` 提取）**。

> 环境说明（2026-11-22 核实）：原始分析环境从未丢失——它以 vendor 目录 `py_lib/` 随项目交付，
> 所有脚本开头 `sys.path.insert(0, .../py_lib)` 都会加载它。系统默认 Python 未装这些包，
> 故直接 `import sklearn` 会失败，但运行项目脚本或设 `PYTHONPATH=py_lib` 即完全可用。
> 已用真实脚本重跑诊断模型、LASSO、网络药理等核心分析，结果与已发表数据**逐位一致（差=0）**。
> 详见 `可复现性验证报告.md`。

**R**：4.4.1（Windows）。由于工作目录含中文路径，需在 C locale 下运行并从 Bash 继承工作目录：

```bash
export LC_ALL=C LANG=C TMPDIR="C:\tmp\rtemp" TMP="C:\tmp\rtemp" TEMP="C:\tmp\rtemp"
"/c/Program Files/R/R-4.4.1/bin/Rscript.exe" scripts/R/reproduce_fig4.R
```

已安装关键包：limma, GEOquery, ggplot2, glmnet, rms, pROC, pheatmap, VennDiagram, corrplot, fgsea。
`clusterProfiler` 系列因需源码编译（本机无 Rtools）未安装，但 ORA 统计已由 Python 完成、GSEA 已由 fgsea 重跑，不影响结论。

## 4. 图表对应关系

- 正文主图：`figures/review/`（Fig1机制 / Fig2流程 / Fig5 AUC / Fig6证据等级 / Fig7覆盖度 / Fig8检查清单），均含 PNG/PDF/SVG。
- 图注：`figures/review/图注草案.md`。
- 案例源图：`figures/Fig*.png`。
- 补充材料 90 项：见 `review_assets/supplementary_index.csv`。

## 5. 引用与数据使用

使用本仓库数据请引用拟投稿的综述论文，并同时引用所涉 GEO 数据集原始文献与关键方法学文献
（完整列表见 `综述正文_第3-4节初稿.md` 附录 A、`review_assets/core_references.ris`）。

## 6. 数据可得性声明模板（可直接用于投稿）

> The expression matrices, analysis scripts, and figures supporting this review are available at https://github.com/xhxccm/pmop-disulfidptosis-review .
> Raw datasets were retrieved from GEO under accession numbers GSE56815, GSE7429, GSE56814, and GSE169396.
> All analyses were conducted in Python 3.13 and R 4.4.1; exact package versions are provided in the repository.
