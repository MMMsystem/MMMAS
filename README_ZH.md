# 孟德尔错配矩阵分析系统 (MMMAS)

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

[English](README.md) | 中文

## 简介

**孟德尔错配矩阵分析系统（Mendelian Mismatch Matrix Analysis System, MMMAS）** 是一款基于布尔孟德尔排除法则的种质库预筛工具。系统通过构建 N×N 错配矩阵，实现种质资源的重复检测、分级评价、遗传兼容性分析等功能。

主要功能包括：

- 自动识别多种 SSR/SNP 基因型数据格式
- 高性能计算孟德尔错配矩阵 M、有效位点矩阵 L、错配率矩阵 MMR
- 自动检测标记类型（SSR / SNP）
- 执行五重数学合规性验证
- 基于 MES（Mendelian Exhaustive Stratification）的分级分析
- 零错配伙伴识别与重复样本检测
- 输出完整的分析报表与可视化数据
- 提供图形界面（GUI）与命令行两种操作模式

---

## 核心概念与术语

### 基础指标

| 缩写 | 英文全称 | 中文名称 | 定义 |
|---|---|---|---|
| MMM | Mendelian Mismatch Matrix | 孟德尔错配矩阵 | N×N 矩阵，元素 M[i,j] 表示样本 i 与 j 的孟德尔不兼容位点数 |
| Mmn | Mendelian Mismatch Number | 孟德尔错配数 | 两个特定样本间的不兼容位点数 |
| Mmin | Mendelian Minimum Mismatch Number | 孟德尔最小错配数 | 某个样本在所有成对比较中的最小错配值 |
| Mmax | Mendelian Maximum Mismatch Number | 孟德尔最大错配数 | 某个样本在所有成对比较中的最大错配值 |
| Mavg | Mendelian Average Mismatch Number | 孟德尔平均错配数 | 某个样本在所有成对比较中的错配平均值 |
| Mzmp | Mendelian Zero-mismatch Partners Number | 孟德尔零错配伙伴数 | 与某样本零错配（Mmn=0）的其他样本数量 |
| Mmmd | Mendelian Mismatch Mode Duplication | 孟德尔错配模式重复 | 零错配对是否为完全重复（1=重复，0=非重复） |
| L | Mendelian Effective Loci Matrix | 孟德尔有效位点矩阵 | N×N 矩阵，记录每对样本间可有效比较的位点数 |
| MMR | Mendelian Mismatch Rate Matrix | 孟德尔错配率矩阵 | N×N 矩阵，记录每对样本间孟德尔错配率（Mmn/L） |

### 分级体系

| 缩写 | 英文全称 | 中文名称 | 定义 |
|---|---|---|---|
| MES | Mendelian Exhaustive Stratification | 孟德尔穷尽分层 | 基于 Mmin 的种质分级体系 |
| Dk | Difference Grade k | 差级 k | Mmin=k 的样本集合 |
| Ck | Cumulative Grade k | 累级 k | Mmin≥k 的样本集合 |
| MGI | Mendelian Grade Gap Index | 孟德尔等级间隙指数 | 空 Dk 等级占总等级的比例 |

---

## 安装与运行环境

### 系统要求

| 项目 | 要求 |
|---|---|
| Python | ≥ 3.8 |
| 操作系统 | Windows / macOS / Linux |
| 内存 | ≥ 4GB（推荐 8GB 以上） |

### 依赖库

```bash
pip install numpy>=1.20.0 pandas>=1.3.0 matplotlib networkx psutil
```

### 安装

解压软件包后直接进入目录运行：

```bash
python MMMAS.py
```

---

## 快速开始

### 方式一：图形界面（推荐初学者）

无参数运行将自动启动 GUI 界面：

```bash
python MMMAS.py
```

操作步骤：

1. 点击“浏览”选择输入 CSV 文件
2. 选择输出目录（默认 `./MMM_Output/`）
3. 选择标记类型（可选，默认自动检测）
4. 设置高级选项（可选）
5. 点击“开始分析”
6. 查看进度和结果

### 方式二：命令行模式

#### 基本用法

```bash
python MMMAS.py -i <输入文件> -o <输出目录> [选项]
```

#### 常用示例

```bash
# 基本用法
python MMMAS.py -i apple_genotype.csv -o ./apple_results/

# 指定 SNP 数据
python MMMAS.py -i snp_data.csv -m SNP -o ./snp_results/

# 大规模数据（>1000 样本）
python MMMAS.py -i large_dataset.csv -o ./large_results/ -cs 500

# 英文界面
python MMMAS.py -i data.csv -o ./results/ -l en
```

### 命令行参数说明

| 参数 | 简写 | 说明 | 默认值 |
|---|---|---|---|
| `--input` | `-i` | 输入基因型数据文件（CSV） | 必填 |
| `--output` | `-o` | 输出目录 | `./MMM_Output/` |
| `--marker-type` | `-m` | 标记类型（SSR/SNP） | 自动检测 |
| `--no-grading` | `-ng` | 关闭自动分级 | False |
| `--no-checkpoint` | `-nc` | 关闭检查点 | False |
| `--no-multithreading` | `-nm` | 关闭多线程 | False |
| `--chunk-size` | `-cs` | 分块大小 | 1000 |
| `--language` | `-l` | 界面语言（zh/en） | zh |
| `--version` | `-v` | 显示版本信息 | - |

---

## 输入数据格式

MMMAS 支持以下数据输入格式，均可自动检测：

### 格式 A：单栏等位基因分离格式（推荐）

```csv
SampleID,Locus1_1,Locus1_2,Locus2_1,Locus2_2,...
Sample1,100,102,150,154,...
Sample2,100,104,150,150,...
```

### 格式 B：单栏基因型格式

```csv
SampleID,Locus1,Locus2,Locus3,...
Sample1,100/102,150/154,200/204,...
Sample2,100/104,150/150,200/202,...
```

### 格式 C：SSR 原始数据格式（双表头）

```csv
SampleID,Locus1,Locus1,Locus2,Locus2,...
,Allele1,Allele2,Allele1,Allele2,...
Sample1,100,102,150,154,...
Sample2,100,104,150,150,...
```

### 格式 D：SNP 转置格式

```csv
SampleID,Sample1,Sample2,Sample3,...
Marker1,AA,AT,TT,...
Marker2,CC,CG,GG,...
```

### SNP 0/1/2 编码格式

```csv
SampleID,SNP1,SNP2,SNP3,...
Sample1,0,1,2
Sample2,1,2,0
```

### 数据要求

| 项目 | 要求 |
|---|---|
| 样本数 | ≥ 2 |
| 位点数 | ≥ 5（推荐 ≥ 11） |
| 缺失数据 | 允许，用 0 或空值表示 |
| 编码方式 | 自动检测 SSR/SNP |

---

## 输出文件说明

运行完成后，输出目录将生成以下文件：

### 核心矩阵文件

| 文件名 | 中文名 | 说明 |
|---|---|---|
| `MMM_Matrix.csv` | 错配矩阵 | N×N 对称矩阵，记录样本 i 与 j 的孟德尔不兼容位点数，对角线为 NA |
| `MMM_Effective_Loci_Matrix.csv` | 有效位点矩阵 | N×N 矩阵，记录每对样本间可有效比较的位点数 |
| `MMM_Mmr_Matrix.csv` | 孟德尔错配率矩阵 | N×N 矩阵，记录每对样本间孟德尔错配率（Mmn/L） |

### 分级文件

| 文件名 | 中文名 | 说明 |
|---|---|---|
| `MMM_basic.csv` | 基础指标表 | 每个样本的核心指标：Mmin、Dk 级别、Mmax、Mavg、Mzmp、Mzmp list. |
| `MMM_Grading_Summary.csv` | 分级汇总表 | 各等级 Dk/Ck 的样本数、占比及成员列表 |

### 零错配与重复文件

| 文件名 | 中文名 | 说明 |
|---|---|---|
| `MMM_Zero_Mismatch_Pairs.csv` | 零错配对 | Mzmp=0 的样本对及 Mmmd 标记 |
| `MMM_Duplicate_Samples.csv` | 重复样本对 | 仅包含 Mmmd=1 的记录 |
| `MMM_Duplicate_Groups.csv` | 重复样本组 | 使用并查集合并的重复样本组 |
| `MMM_only_data.csv` | 去重后数据 | 自动去除重复样本后的基因型数据 |

### 频率统计文件

| 文件名 | 中文名 | 说明 |
|---|---|---|
| `MMM_Mnp_Mismatch_Pairs.csv` | 低错配品种对 | 错配数在 0~2 之间的品种对列表 |
| `MMM_frequency.csv` | 频率分布 | Mmn、Mmin、Mzmp、Mavg 的分布频率 |

### 报告文件

| 文件名 | 中文名 | 说明 |
|---|---|---|
| `MMM_Analytical_Report.csv` | 分析报告 | 包含运行参数、验证结果、分级概要的表格报告 |
| `MMM_Build_Report.txt` | 构建报告 | 包含运行时间、各阶段耗时、输出文件列表、文件说明的文本报告 |
| `MMM_MGI.txt` | 孟德尔等级间隙指数 | MGI、Max Mmin、Total Samples 等信息 |

---

## 核心模块架构

```
MMMAS.py
├── MMMConfig              # 配置管理类
├── MMMResults             # 结果容器类
├── DataStructureDetector  # 数据格式自动检测
├── MarkerTypeDetector     # 标记类型自动检测
├── GenotypeProcessor      # 基因型数组构建
├── MismatchCalculator     # 错配计算（标准/向量化）
├── ParallelMMMCalculator  # 并行错配计算类
├── MMMValidator           # 矩阵合规性验证类
├── DkClassifier           # MES 分级类
├── GradingAnalyzer        # 分级分析类
├── ReportGenerator        # 报告生成类
└── MMMWorkflow            # 统一工作流类
```

工作流程：

```
输入 CSV → 格式检测 → 基因型数组构建 → 错配矩阵计算 → 矩阵验证 → Dk 分级 → 零错配检测 → 报告生成
```

---

## 算法原理

### 孟德尔不兼容判定

对于两位点基因型 `(a,b)` 和 `(c,d)`，孟德尔不兼容判定规则：

- 如果两个基因型共享至少一个等位基因，则兼容（M=0）
- 如果两个基因型无共享等位基因，则不兼容（M=1）

形式化表达：

- 兼容条件：`{a,b} ∩ {c,d} ≠ ∅`
- 不兼容条件：`{a,b} ∩ {c,d} = ∅`

### 错配数计算

对于两个样本 i 和 j，错配数 M[i,j] 为所有位点上孟德尔不兼容判定的累加：

```
M[i,j] = Σ(Incompatible(G[i,l], G[j,l]))  对于所有有效位点 l
```

其中 `Incompatible()` 在不兼容时返回 1，兼容时返回 0。

### MES 分级算法

- 计算每个样本的 Mmin（行最小值，排除对角线）
- `Dk = {样本 | Mmin = k}`
- `Ck = {样本 | Mmin ≥ k}`
- `MGI = 空 Dk 等级数 / (最大 Mmin + 1)`

### 重复检测算法

1. 找出所有 Mmn=0 的样本对
2. 对每对零错配样本，比较它们与所有第三方样本的错配轮廓
3. 如果两样本的错配轮廓完全相同 → Mmmd=1（重复）
4. 使用并查集将重复对合并为重复组

---

## 五重数学验证

构建的错配矩阵自动通过以下合规性验证：

| 验证项 | 说明 |
|---|---|
| 对称性 | 矩阵满足 `M[i,j] = M[j,i]` |
| 非负性 | 所有错配数 ≥ 0 |
| 有界性 | 所有错配数 ≤ 总位点数 |
| 整数性 | 所有错配数为整数 |
| 对角线规则 | 对角线元素为缺失值（NaN） |

---

## 性能优化

### 计算策略选择

| 样本数 | 策略 | 时间复杂度 |
|---|---|---|
| ≤ 100 | 标准循环 | O(N² × L) |
| ≤ 1,000 | 向量化计算 | O(N² × L) |
| > 1,000 | 并行分块计算 | O(N² × L / P) |

### 性能参考

| 样本数 | 位点数 | 预估时间 | 内存 |
|---|---|---|---|
| 100 | 15 | < 1 秒 | ~50MB |
| 500 | 15 | ~5 秒 | ~200MB |
| 1,000 | 15 | ~30 秒 | ~800MB |
| 5,000 | 15 | ~5 分钟 | ~4GB |

---

## 常见问题

### Q1: 如何处理大量样本（>5000）？

建议开启多线程并调整分块大小：

```bash
python MMMAS.py -i large.csv -o ./output/ -cs 500
```

### Q2: Mzmp 值很高说明什么？

Mzmp（零错配伙伴数）高说明该样本在库中有多遗传相同的品种，可能为：

- 广泛栽培的主栽品种
- 历史育种核心亲本
- 同物异名品种

### Q3: Mmmd=0 和 Mmmd=1 的区别？

- **Mmmd=1**：两个样本与所有其他样本的错配模式完全相同，提示为重复样本
- **Mmmd=0**：两个样本仅彼此零错配，但与第三方样本的错配模式不同，可能是极近亲缘关系

### Q4: MGI 值如何解读？

| MGI 范围 | 解读 |
|---|---|
| MGI = 0 | 等级连续，无间隙（理想状态） |
| MGI < 0.1 | 等级基本连续 |
| MGI 0.1~0.3 | 存在明显间隙 |
| MGI > 0.3 | 等级严重不连续，可能需要检查数据质量 |

### Q5: Dk 分级的生物学意义？

| Dk 级别 | 生物学意义 |
|---|---|
| Dk=0 | 核心种质，库内无更近亲缘 |
| Dk=1 | 至少有一个零错配伙伴，存在近缘关系 |
| Dk≥2 | 遗传距离逐步增大，亲缘关系疏远 |

---

## 版本信息

- **当前版本**：v1.0.0
- **发布日期**：2026-07-19
- **作者**： 陈启亮、陈刘秀等

---

## 许可证

本项目采用 MIT 许可证开源，详见 [LICENSE](LICENSE) 文件。

---

## 引用

如果您在研究中使用了本系统，请引用：

> MMM Assistant. (2026). Mendelian Mismatch Matrix Analysis System (MMM v1.0).

---

## 联系方式

| 项目 | 信息 |
|---|---|
| 电子邮箱 | 26451851@qq.com|
| 单位/机构 | 湖北省农业科学院果树茶叶研究所|

如有问题或建议，欢迎通过上述方式联系。
