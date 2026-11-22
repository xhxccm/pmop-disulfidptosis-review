# 第七阶段：实验验证方案（线下执行计划）

> 说明：本阶段为湿实验，无法在计算环境执行。以下为基于前面生物信息学结果设计的实验方案，
> 供线下开展。核心 Hub Genes：ADRB1、FLT3、MAP3K1、MEX3D、OSTF1、DHX35 等 18 个；
> 机制核心（DRG-DEG 交集）：FLNC、PGRMC2；蒙药蓝刺头重点成分：槲皮素、芹菜素、咖啡酸、绿原酸。

---

## 7.1 细胞实验：成骨细胞双硫死亡模型

**目的**：在成骨细胞（MC3T3-E1 或 hFOB1.19）中建立双硫死亡模型，验证蓝刺头成分的保护作用。

**材料**
- 细胞：MC3T3-E1（小鼠成骨前体）/ hFOB1.19（人成骨）
- 诱导剂：葡萄糖转运抑制（2-DG）联合胱氨酸缺乏培养基，或 SLC7A11 抑制剂（如 Erastin 在特定条件下诱导双硫死亡）
- 药物：槲皮素、芹菜素（Sigma），浓度梯度 5/10/20 μM
- 载体：青链霉素、DMEM/α-MEM、胎牛血清

**分组**
1. 对照组（正常培养基）
2. 模型组（双硫死亡诱导）
3. 模型 + 蓝刺头成分低/中/高剂量组
4. 阳性对照（NAC 或 DTT，抗氧化剂）

**指标**
- 细胞活力：CCK-8 / LDH 释放
- 双硫死亡形态：F-actin 免疫荧光（鬼笔环肽）观察肌动蛋白网络收缩/塌陷
- 氧化应激：ROS（DCFH-DA）、GSH/GSSG 比值、NADP+/NADPH
- 双硫死亡标志：SLC7A11、ACTB（β-actin）聚集状态、WRC 复合体（RAC1、NCKAP1、CYFIP1）

---

## 7.2 Western Blot：SLC7A11、ACTB、GSH 相关蛋白

**目的**：验证双硫死亡通路关键蛋白及 Hub Gene 蛋白表达。

**检测蛋白**
- 双硫死亡核心：SLC7A11、ACTB、RAC1、NCKAP1、G6PD、TXNRD1
- Hub 蛋白：ADRB1、MAP3K1、PGRMC2、FLNC
- 内参：GAPDH / β-tubulin

**方法要点**
- RIPA 裂解 + BCA 定量；10–12% SDS-PAGE；转 PVDF 膜
- 一抗 4°C 过夜，二抗室温 1 h；ECL 显影，ImageJ 灰度定量
- 每组 3 个生物学重复，单因素 ANOVA + Tukey 事后检验

> 注：GSH 本身为小分子，用 GSH/GSSG 检测试剂盒（酶循环法）定量，不做 WB。

---

## 7.3 RT-qPCR：Hub Genes mRNA 表达验证

**目的**：在细胞模型及（如有）临床样本中验证 Hub Genes 转录水平。

**基因**：ADRB1、FLT3、MAP3K1、MEX3D、OSTF1、DHX35、FLNC、PGRMC2（引物用 Primer3 设计，跨外显子）
**内参**：GAPDH、ACTB
**方法**：Trizol 提 RNA → 逆转录 cDNA → SYBR Green qPCR（ABI 7500），2^-ΔΔCt 法
**重复**：每样本 3 技术重复，≥3 生物学重复

---

## 7.4 动物实验：去卵巢（OVX）大鼠模型 + 蒙药干预

**目的**：在体验证蓝刺头对绝经后骨质疏松的保护作用及双硫死亡机制。

**动物**：SD 大鼠，雌性，12 周龄，SPF 级
**分组**（n=8/组）
1. Sham（假手术）
2. OVX（去卵巢）
3. OVX + 蓝刺头提取物低剂量
4. OVX + 蓝刺头提取物高剂量
5. OVX + 阳性药（如雌二醇或阿仑膦酸钠）

**流程**
- 双侧卵巢切除建立 PMOP 模型，术后 4 周开始灌胃给药，连续 12 周
- 处死前取材：股骨/胫骨（骨密度、Micro-CT、三点弯曲）、血清（PINP、CTX、ROS、GSH）
- 骨组织：HE/TRAP 染色（破骨细胞）、免疫组化（ADRB1、PGRMC2、SLC7A11）、Western Blot / qPCR

**主要终点**
- Micro-CT：BV/TV、Tb.N、Tb.Th、Tb.Sp
- 双硫死亡相关：骨组织 F-actin 形态、SLC7A11/ACTB 表达、GSH 水平

---

## 关键试剂与指标清单

| 类别 | 项目 | 用途 |
|---|---|---|
| 诱导剂 | 2-DG、胱氨酸缺乏培养基、Erastin | 双硫死亡建模 |
| 药物 | 槲皮素、芹菜素、咖啡酸、绿原酸 | 蓝刺头单体干预 |
| 检测 | CCK-8、LDH、DCFH-DA、GSH/GSSG 试剂盒 | 活力/氧化应激 |
| 抗体 | SLC7A11、ACTB、RAC1、ADRB1、PGRMC2、FLNC | WB/IHC |
| 分子 | qPCR 引物（8 Hub + 内参） | mRNA 验证 |
| 动物 | SD 大鼠、雌二醇/阿仑膦酸钠 | OVX 模型与阳性对照 |
| 影像 | Micro-CT | 骨微结构 |

---

## 与计算结果的衔接

- 若细胞/动物实验中双硫死亡指标（F-actin 塌陷、GSH 耗竭）被蓝刺头成分逆转，且 Hub 蛋白表达随之改变，
  则支持"蒙药蓝刺头通过调控双硫死亡干预 PMOP"的核心假说。
- 分子对接显示槲皮素/芹菜素与 MAP3K1、FLT3、ADRB1、PGRMC2 结合能 −7.3 ~ −8.3 kcal/mol，
  可在 WB/qPCR 中优先验证这些靶点。
