# Methodological review companion repository

Companion code, data, figures and manuscript drafts for:/n*Pipelines, Pitfalls, and External Validation in Multi-omics Biomarker Discovery: A Worked Example from Disulfidptosis-Related Genes in Postmenopausal Osteoporosis*.

## Contents
- `scripts/python/`, `scripts/R/` — all analysis and plotting scripts
- `data/` — expression matrices and all result tables (CSVs)
- `figures/main/` — Figures 1-8 (PNG + SVG)
- `manuscript/` — review framework, full draft sections, figure legends
- `full_requirements_exact.txt` — 64 packages with exact versions
- `environment_R_session.txt` — R 4.4.1 package versions
- `可复现性验证报告.md` — 13-item zero-difference re-run verification

## Reproduce
Set `PYTHONPATH` to a `py_lib` environment (see 可复现性验证报告.md) or install
the exact versions from `full_requirements_exact.txt`, then run any script in `scripts/`.

## Data availability
Raw data from GEO: GSE56815, GSE7429, GSE56814, GSE169396.
