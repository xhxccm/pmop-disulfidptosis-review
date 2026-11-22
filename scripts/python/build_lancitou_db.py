"""
蒙药蓝刺头（Echinops latifolius / Echinops sphaerocephalus）活性成分-靶点数据库构建

证据分级：
  A - 实验证据：蓝刺头提取物/含药血清在动物或细胞实验中直接检测到靶点表达变化
  B - 策展证据：CTD (Comparative Toxicogenomics Database) curated chemical-gene interactions
  C - 预测证据：蓝刺头文献中的网络药理学/分子对接预测（未经该药材实验验证）

输出：
  mongolian_lancitou_compound_target.csv  成分-靶点主表
  mongolian_lancitou_compounds.csv        成分信息表
  mongolian_lancitou_targets.csv          靶点信息表
  mongolian_lancitou_ct.db                SQLite 数据库
  build_lancitou_db.py                   可复现脚本
"""

import io
import pathlib
import sqlite3
import time

import pandas as pd
import requests

BASE = pathlib.Path(__file__).resolve().parent

COMPOUNDS = {
    "quercetin": {"cn": "槲皮素", "type": "单体"},
    "apigenin": {"cn": "芹菜素（4',5,7-三羟基黄酮）", "type": "单体"},
    "caffeic acid": {"cn": "咖啡酸", "type": "单体"},
    "chlorogenic acid": {"cn": "绿原酸", "type": "单体"},
    "palmitic acid": {"cn": "正十六烷酸", "type": "单体"},
}

HERB_CN = "蒙药蓝刺头"
HERB_EN = "Echinops latifolius / Echinops sphaerocephalus (Mongolian medicine)"

# 蓝刺头研究中有明确文献来源的实验证据（A 级）
STUDY_A = [
    # 抗骨质疏松研究：去卵巢 + D-半乳糖大鼠，ELISA 检测血清靶点表达
    ("蓝刺头总酚酸/水提物", "Lancitou total phenolic acids / aqueous extract", "混合物", "HSD11B1",
     "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "去卵巢及D-半乳糖骨质疏松大鼠，低/中/高剂量均下调血清HSD11B1", "PMOP rat model"),
    ("蓝刺头总酚酸/水提物", "Lancitou total phenolic acids / aqueous extract", "混合物", "CYP19A1",
     "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "去卵巢及D-半乳糖骨质疏松大鼠，低/中/高剂量均下调血清CYP19A1", "PMOP rat model"),
    ("蓝刺头总酚酸/水提物", "Lancitou total phenolic acids / aqueous extract", "混合物", "ALOX15",
     "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "去卵巢及D-半乳糖骨质疏松大鼠，低/中/高剂量均下调血清ALOX15", "PMOP rat model"),
    # 促血管生成研究：HUVEC qPCR
    ("蓝刺头乙醇提取物", "Lancitou ethanol extract", "混合物", "HIF1A",
     "方源等. 中国组织工程研究 2025;29(35):7519-7528",
     "HUVEC 10 μg/mL 干预，qPCR 显示 HIF-1α mRNA 上调", "HUVEC angiogenesis"),
    ("蓝刺头乙醇提取物", "Lancitou ethanol extract", "混合物", "VEGFA",
     "方源等. 中国组织工程研究 2025;29(35):7519-7528",
     "HUVEC 10 μg/mL 干预，qPCR 显示 VEGFA mRNA 上调", "HUVEC angiogenesis"),
    ("蓝刺头乙醇提取物", "Lancitou ethanol extract", "混合物", "KDR",
     "方源等. 中国组织工程研究 2025;29(35):7519-7528",
     "HUVEC 10 μg/mL 干预，qPCR 显示 KDR mRNA 上调", "HUVEC angiogenesis"),
    # PMOP BMSCs 研究：RT-PCR / Western blot
    ("蓝刺头含药血清", "Lancitou medicated serum", "混合物", "WNT3A",
     "杨嘉琳. 内蒙古医科大学硕士论文 2023（Wnt/β-catenin 通路）",
     "BMSCs 含药血清干预，Wnt3a mRNA 与蛋白升高；DKK-1 可抑制", "PMOP BMSCs"),
    ("蓝刺头含药血清", "Lancitou medicated serum", "混合物", "CTNNB1",
     "杨嘉琳. 内蒙古医科大学硕士论文 2023（Wnt/β-catenin 通路）",
     "BMSCs 含药血清干预，β-catenin mRNA 与蛋白升高；DKK-1 可抑制", "PMOP BMSCs"),
    ("蓝刺头含药血清", "Lancitou medicated serum", "混合物", "RUNX2",
     "杨嘉琳. 内蒙古医科大学硕士论文 2023（Wnt/β-catenin 通路）",
     "BMSCs 含药血清干预，Runx2 mRNA 与蛋白升高；DKK-1 可抑制", "PMOP BMSCs"),
]

# 蓝刺头文献中的网络药理学预测（C 级）
STUDY_C = [
    ("槲皮素", "quercetin", "HSD11B1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("槲皮素", "quercetin", "CYP19A1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("槲皮素", "quercetin", "ALOX15", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("芹菜素", "apigenin", "HSD11B1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("芹菜素", "apigenin", "CYP19A1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("芹菜素", "apigenin", "ALOX15", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("咖啡酸", "caffeic acid", "HSD11B1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("咖啡酸", "caffeic acid", "CYP19A1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("咖啡酸", "caffeic acid", "ALOX15", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("正十六烷酸", "palmitic acid", "HSD11B1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("正十六烷酸", "palmitic acid", "CYP19A1", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
    ("正十六烷酸", "palmitic acid", "ALOX15", "基于网络药理学研究蒙药蓝刺头抗骨质疏松的作用机制",
     "Swiss/SuperPred 预测 + 网络拓扑筛选的关键成分-靶点"),
]

ANGIO_CORE_GENES = ["TGFB1", "TNF", "IL6", "STAT3", "CTNNB1", "IL1B", "AKT1", "HIF1A"]
for cn, en, _, _ in [
    ("槲皮素", "quercetin", "", ""), ("芹菜素", "apigenin", "", ""),
    ("咖啡酸", "caffeic acid", "", ""), ("绿原酸", "chlorogenic acid", "", ""),
]:
    for g in ANGIO_CORE_GENES:
        STUDY_C.append((
            cn, en, g, "方源等. 中国组织工程研究 2025;29(35):7519-7528",
            "网络药理学核心活性成分-核心靶点（促血管生成）",
        ))

PMOP_CORE_GENES = ["ALB", "AKT1", "CTNNB1", "EGFR"]
for g in PMOP_CORE_GENES:
    STUDY_C.append((
        "蓝刺头主要成分（网络药理学筛选）", "Lancitou major constituents (screened)", g,
        "杨嘉琳. 内蒙古医科大学硕士论文 2023（Wnt/β-catenin 通路）",
        "PharmMapper + 疾病数据库交集(PPI)筛选的核心靶点",
    ))


def fetch_pubchem(name: str) -> dict:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,MolecularFormula/JSON"
    try:
        r = requests.get(url, timeout=(15, 45))
        if r.status_code != 200:
            return {}
        props = r.json()["PropertyTable"]["Properties"][0]
        return {
            "pubchem_cid": props.get("CID"),
            "molecular_formula": props.get("MolecularFormula"),
            "smiles": props.get("ConnectivitySMILES") or props.get("CanonicalSMILES"),
        }
    except Exception:
        return {}


def fetch_ctd(name: str) -> pd.DataFrame:
    url = "https://ctdbase.org/tools/batchQuery.go"
    params = {"inputType": "chem", "inputTerms": name, "report": "genes_curated", "format": "tsv"}
    for attempt in range(3):
        try:
            r = requests.get(url, params=params, timeout=(20, 120))
            r.raise_for_status()
            df = pd.read_csv(io.StringIO(r.text), sep="\t")
            df = df[df["Organism"] == "Homo sapiens"].copy()
            df["input_name"] = name
            return df
        except Exception as e:
            if attempt == 2:
                print(f"  CTD 查询失败 {name}: {e}")
                return pd.DataFrame()
            time.sleep(2)
    return pd.DataFrame()


def main():
    print("[1/5] 获取化合物元数据 (PubChem) ...")
    compound_rows = []
    meta = {}
    for en, info in COMPOUNDS.items():
        m = fetch_pubchem(en)
        meta[en] = m
        compound_rows.append({
            "name_cn": info["cn"], "name_en": en, "compound_type": info["type"],
            "pubchem_cid": m.get("pubchem_cid"), "molecular_formula": m.get("molecular_formula"),
            "smiles": m.get("smiles"), "ctd_chemical_id": None, "cas_rn": None,
            "evidence_note": "蓝刺头文献中明确提及的活性成分",
        })
        print(f"      {en}: CID={m.get('pubchem_cid')}")
    for cn, en, ctype in [
        ("蓝刺头总酚酸/水提物", "Lancitou total phenolic acids / aqueous extract", "混合物"),
        ("蓝刺头乙醇提取物", "Lancitou ethanol extract", "混合物"),
        ("蓝刺头含药血清", "Lancitou medicated serum", "混合物"),
        ("蓝刺头主要成分（网络药理学筛选）", "Lancitou major constituents (screened)", "成分集合"),
    ]:
        compound_rows.append({
            "name_cn": cn, "name_en": en, "compound_type": ctype,
            "pubchem_cid": None, "molecular_formula": None, "smiles": None,
            "ctd_chemical_id": None, "cas_rn": None,
            "evidence_note": "来自蓝刺头整体/成分集合的文献证据，非单一化合物",
        })
    compounds_df = pd.DataFrame(compound_rows)

    print("[2/5] 获取 CTD curated 成分-基因互作 (Level B) ...")
    ctd_frames = []
    for en in COMPOUNDS:
        df = fetch_ctd(en)
        if not df.empty:
            ctd_frames.append(df)
            print(f"      {en}: {len(df)} 条人类互作")
    ctd = pd.concat(ctd_frames, ignore_index=True) if ctd_frames else pd.DataFrame()

    # 回填化合物表的 CTD ID / CAS
    if not ctd.empty:
        ctd_map = ctd.groupby("input_name").agg(
            ctd_chemical_id=("ChemicalId", "first"), cas_rn=("CasRN", "first")
        ).to_dict("index")
        for en, vals in ctd_map.items():
            mask = compounds_df["name_en"] == en
            compounds_df.loc[mask, "ctd_chemical_id"] = vals["ctd_chemical_id"]
            compounds_df.loc[mask, "cas_rn"] = vals["cas_rn"]

    rows = []
    if not ctd.empty:
        for _, r in ctd.iterrows():
            en = r["input_name"]
            rows.append({
                "herb_cn": HERB_CN, "herb_en": HERB_EN,
                "compound_cn": COMPOUNDS[en]["cn"], "compound_en": en,
                "target_symbol": r["GeneSymbol"], "target_entrez_id": r["GeneId"],
                "evidence_level": "B",
                "evidence_type": "CTD curated chemical-gene interaction",
                "source": "CTD (https://ctdbase.org/tools/batchQuery.go, genes_curated)",
                "source_detail": f"chemical={r['ChemicalName']}, CTD ID={r['ChemicalId']}, CAS={r['CasRN']}",
                "pubmed_ids": r.get("PubMedIds"), "context": "Homo sapiens curated literature",
                "note": "非蓝刺头特异；需结合药材成分存在性判断",
            })

    print("[3/5] 加入蓝刺头研究实验证据 (Level A) ...")
    for cn, en, ctype, gene, source, detail, context in STUDY_A:
        rows.append({
            "herb_cn": HERB_CN, "herb_en": HERB_EN,
            "compound_cn": cn, "compound_en": en,
            "target_symbol": gene, "target_entrez_id": None,
            "evidence_level": "A",
            "evidence_type": "experimental (in vivo / in vitro expression change)",
            "source": source, "source_detail": detail,
            "pubmed_ids": None, "context": context,
            "note": "蓝刺头药材/提取物/含药血清层面的实验证据，未拆到单一成分",
        })

    print("[4/5] 加入蓝刺头网络药理学预测 (Level C) ...")
    for cn, en, gene, source, detail in STUDY_C:
        rows.append({
            "herb_cn": HERB_CN, "herb_en": HERB_EN,
            "compound_cn": cn, "compound_en": en,
            "target_symbol": gene, "target_entrez_id": None,
            "evidence_level": "C",
            "evidence_type": "network pharmacology prediction / molecular docking",
            "source": source, "source_detail": detail,
            "pubmed_ids": None, "context": "predicted, not experimentally validated for this herb",
            "note": "预测证据，优先级低于 A/B",
        })

    cpt_df = pd.DataFrame(rows)
    # 仅去除完全同源的重复记录（同一成分-靶点-来源）；不同证据等级/来源均保留，
    # 以便区分「蓝刺头特异性证据」与「CTD 通用策展证据」。
    cpt_df = cpt_df.drop_duplicates(
        ["compound_en", "target_symbol", "source"], keep="first"
    ).reset_index(drop=True)

    # 用 CTD 映射回填靶点 Entrez ID（若有）
    if not ctd.empty:
        gene2entrez = ctd.drop_duplicates("GeneSymbol").set_index("GeneSymbol")["GeneId"].to_dict()
        cpt_df["target_entrez_id"] = cpt_df.apply(
            lambda r: r["target_entrez_id"] or gene2entrez.get(r["target_symbol"]), axis=1
        )

    # 与 GSE56814 表达矩阵交叉标记
    expr_path = BASE / "GSE56814_exprs_matrix.csv"
    if expr_path.exists():
        expr_genes = set(pd.read_csv(expr_path, index_col=0).index.astype(str))
        cpt_df["in_GSE56814_matrix"] = cpt_df["target_symbol"].isin(expr_genes)
    else:
        cpt_df["in_GSE56814_matrix"] = None

    # 靶点表
    target_rows = []
    if not ctd.empty:
        for gene, grp in ctd.groupby("GeneSymbol"):
            target_rows.append({
                "target_symbol": gene,
                "target_entrez_id": grp["GeneId"].iloc[0],
                "organism": "Homo sapiens",
                "n_ctd_compounds": len(grp),
            })
    targets_df = pd.DataFrame(target_rows)
    if not targets_df.empty:
        targets_df = targets_df.merge(
            cpt_df[["target_symbol", "in_GSE56814_matrix"]].drop_duplicates("target_symbol"),
            on="target_symbol", how="left",
        )

    print("[5/5] 写出文件 ...")
    cpt_df.to_csv(BASE / "mongolian_lancitou_compound_target.csv", index=False)
    compounds_df.to_csv(BASE / "mongolian_lancitou_compounds.csv", index=False)
    targets_df.to_csv(BASE / "mongolian_lancitou_targets.csv", index=False)

    db_path = BASE / "mongolian_lancitou_ct.db"
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(db_path)
    compounds_df.to_sql("compounds", con, index=False)
    cpt_df.to_sql("compound_target", con, index=False)
    targets_df.to_sql("targets", con, index=False)
    con.execute("CREATE INDEX idx_ct_compound ON compound_target(compound_en)")
    con.execute("CREATE INDEX idx_ct_target ON compound_target(target_symbol)")
    con.execute("CREATE INDEX idx_ct_level ON compound_target(evidence_level)")
    con.commit()
    con.close()

    print("完成。")
    print(f"  成分数: {len(compounds_df)}")
    print(f"  成分-靶点记录: {len(cpt_df)}")
    if "in_GSE56814_matrix" in cpt_df:
        print(f"  其中在 GSE56814 表达矩阵中可检测: {int(cpt_df['in_GSE56814_matrix'].sum())}")
    print(f"  证据等级分布:\n{cpt_df['evidence_level'].value_counts().to_string()}")


if __name__ == "__main__":
    main()
