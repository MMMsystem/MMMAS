# 孟德尔错配矩阵分析系统 (MMM)

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## 简介

**孟德尔错配矩阵分析系统 (Mendelian Mismatch Matrix Analysis System, MMM)** 是一款用于植物品种遗传关系分析的专业工具。该系统基于孟德尔遗传定律，通过计算样本间的错配矩阵（Mismatch Matrix），实现对品种间遗传差异的定量评估与分级分类。

本系统整合了**矩阵构建**与**分级分析**两大核心模块，支持 SSR 和 SNP 两种主流分子标记类型，能够自动识别多种基因型数据格式，并提供图形用户界面（GUI）与命令行两种操作模式。

---

## 核心功能

### 1. 自动格式识别与矩阵构建
- 自动识别多种基因型数据格式：
  - `original`：斜杠分隔格式（如 `120/130`）
  - `reformatted`：分列格式（如 `_1`/`_2`、`.1`/`.2`、`_a`/`_b` 后缀）
  - `SNP`：0/1/2 编码格式
  - `transposed_base_snp`：转置碱基型格式（如 `CC`、`TC`、`TT`）
- 自动检测标记类型（SSR / SNP）
- 支持 CSV 文件多种编码自动识别（UTF-8、GBK、Latin-1 等）

### 2. 高性能错配矩阵计算
- **标准算法**：三重循环，适用于 ≤100 个样本的小规模数据
- **向量化算法**：基于 NumPy 矩阵运算，适用于 ≤1000 个样本的中等规模数据
- **并行计算**：多线程分块处理，适用于 >1000 个样本的大规模数据

### 3. 六重数学验证
构建的错配矩阵自动通过以下合规性验证：
| 验证项 | 说明 |
|--------|------|
| 对称性 | 矩阵满足 `M[i,j] = M[j,i]` |
| 非负性 | 所有错配数 ≥ 0 |
| 有界性 | 所有错配数 ≤ 总位点数 |
| 整数性 | 所有错配数为整数 |
| 对角线规则 | 对角线元素为缺失值（NaN） |
| 三角不等式 | 满足 `M[i,j] + M[j,k] ≥ M[i,k]` |

### 4. 分级分析体系
- **基础指标计算**：Mmin（最小错配数）、Mmax（最大错配数）、Mavg（平均错配数）、Mzmp（零错配伙伴数）
- **差级分类 (Dk)**：按 Mmin 值将样本划分为不同差异等级
- **累级分类 (Gk)**：按 Mmin ≥ k 的累积关系构建层级体系
- **间隙指数 (MGI)**：评估分级体系的连续性
- **零错配对检测**：识别孟德尔零错配样本对
- **重复样本检测 (Mmmd)**：基于并查集算法检测完全重复的样本组

### 5. D0 网络连通分析
- 针对 Mmin=0 的样本构建连通网络
- 识别连通分量、孤立样本及网络拓扑特征

### 6. 图形用户界面 (GUI)
- 基于 Tkinter 的友好操作界面
- **矩阵构建标签页**：可视化选择输入文件、设置参数、实时查看日志
- **分级分析标签页**：提供样本筛选、样本对比、样本搜索、分级体系查看、重复样本检测、D0 网络连通分析等功能

---

## 安装依赖

```bash
pip install numpy pandas matplotlib networkx
```

### 系统要求
- Python ≥ 3.8
- 操作系统：Windows / macOS / Linux

---

## 快速开始

### 方式一：图形界面（推荐初学者）

```bash
python MMMv1.py
```

无参数运行将自动启动 GUI 界面：
1. 选择**输入 CSV 文件**
2. 选择**输出目录**
3. 点击【开始构建】

### 方式二：命令行模式

#### 构建矩阵并自动分级
```bash
python MMMv1.py -i input.csv -o output_dir
```

#### 仅对已有矩阵进行分级分析
```bash
python MMMv1.py --grading -m MMM_Matrix.csv -o grading_output
```

#### 指定 SSR 标记类型
```bash
python MMMv1.py -i input.csv -o output_dir -t SSR
```

#### 禁用多线程（调试时使用）
```bash
python MMMv1.py -i input.csv -o output_dir --no-threads
```

### 命令行参数说明

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入原始基因型 CSV 文件路径 |
| `-o, --output` | 输出目录路径 |
| `-m, --matrix` | 已有 MMM 矩阵文件路径（分级模式） |
| `-t, --type` | 标记类型：`SSR` 或 `SNP`（默认自动检测） |
| `--no-grading` | 构建矩阵后不自动进行分级分析 |
| `--grading` | 仅执行分级分析模式 |
| `--no-threads` | 禁用多线程 |
| `-v, --version` | 显示版本信息 |

---

## 输入数据格式示例

### 格式 1：Original（斜杠分隔）
```csv
SampleID,LOCUS1,LOCUS2,LOCUS3
Sample001,120/130,140/150,160/170
Sample002,120/120,140/140,160/160
```

### 格式 2：Reformatted（分列）
```csv
SampleID,LOCUS1_1,LOCUS1_2,LOCUS2_1,LOCUS2_2
Sample001,120,130,140,150
Sample002,120,120,140,140
```

### 格式 3：SNP（0/1/2 编码）
```csv
SampleID,SNP1,SNP2,SNP3
Sample001,0,1,2
Sample002,0,0,1
```

---

## 输出文件说明

运行完成后，输出目录将生成以下文件：

| 文件名 | 说明 |
|--------|------|
| `MMM_Matrix.csv` | 错配矩阵 M |
| `MMM_Effective_Loci_Matrix.csv` | 有效位点矩阵 L |
| `MMM_Analytical_Report.csv` | 完整分析报告 |
| `MMM_Dk_Classification.csv` | Dk 分级结果表 |
| `MMM_basic.csv` | 基础指标表（Mmin/Mmax/Mavg/Mzmp） |
| `MMM_Grading_Summary.csv` | 分级体系汇总表 |
| `MMM_MGI.txt` | 间隙指数信息 |
| `MMM_Zero_Mismatch_Pairs.csv` | 零错配对检测结果 |
| `MMM_Duplicate_Samples.csv` | 重复样本对列表 |
| `MMM_Duplicate_Groups.csv` | 重复样本组列表 |
| `MMM_Build_Report.txt` | 构建摘要报告（文本格式） |

---

## 算法原理

### 孟德尔错配判定
对于两个样本的基因型 `g1 = (a1, a2)` 和 `g2 = (b1, b2)`，当它们**没有任何共享等位基因**时，判定为孟德尔错配：

```
is_mismatch(g1, g2) = true  当且仅当  {a1, a2} ∩ {b1, b2} = ∅
```

### 错配矩阵定义
对于 N 个样本、L 个位点的数据集，错配矩阵 M 是一个 N×N 的对称矩阵：

```
M[i,j] = 样本 i 与样本 j 之间存在孟德尔错配的位点数
```

对角线元素 `M[i,i]` 设为缺失值（NaN）。

### Dk 差级分类
样本 `s` 的差级 `Dk` 定义为其与其他所有样本的最小错配数：

```
Dk(s) = min{ M[s, j] | j ≠ s }
```

差级越小，表示该样本在群体中具有越近的遗传关系。

---

## 项目结构

```
.
├── MMMv1.py          # 主程序（中文版）
├── MMMv1_en.py       # 主程序（英文版）
├── README.md         # 中文说明文档
├── README_EN.md      # 英文说明文档
└── ...
```

---

## 版本信息

- **当前版本**：v1.0.0
- **发布日期**：2026-04-19
- **作者**：MMM Assistant

---

## 许可证

本项目采用 MIT 许可证开源，详见 [LICENSE](LICENSE) 文件。

---

## 引用与致谢

如果您在研究中使用了本系统，请引用：

> MMM Assistant. (2026). Mendelian Mismatch Matrix Analysis System (MMM v1.0). 

---

## 联系我们

如有问题或建议，欢迎通过 GitHub Issues 提交反馈。
