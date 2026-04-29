# Mendelian Mismatch Matrix Analysis System (MMM)

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Introduction

The **Mendelian Mismatch Matrix Analysis System (MMM)** is a professional tool for analyzing genetic relationships among plant varieties. Based on Mendelian inheritance laws, this system calculates the mismatch matrix between samples to quantitatively evaluate and classify genetic differences.

The system integrates two core modules: **Matrix Construction** and **Grading Analysis**, supports both SSR and SNP molecular marker types, automatically recognizes multiple genotype data formats, and provides both Graphical User Interface (GUI) and command-line operation modes.

---

## Core Features

### 1. Automatic Format Recognition and Matrix Construction
- Automatically recognizes multiple genotype data formats:
  - `original`: Slash-delimited format (e.g., `120/130`)
  - `reformatted`: Split-column format (e.g., `_1`/`_2`, `.1`/`.2`, `_a`/`_b` suffixes)
  - `SNP`: 0/1/2 coded format
  - `transposed_base_snp`: Transposed base format (e.g., `CC`, `TC`, `TT`)
- Automatic marker type detection (SSR / SNP)
- Automatic CSV encoding recognition (UTF-8, GBK, Latin-1, etc.)

### 2. High-Performance Mismatch Matrix Calculation
- **Standard Algorithm**: Triple-loop implementation, suitable for small datasets (≤100 samples)
- **Vectorized Algorithm**: NumPy-based matrix operations, suitable for medium datasets (≤1000 samples)
- **Parallel Computing**: Multi-threaded block processing, suitable for large datasets (>1000 samples)

### 3. Six-Fold Mathematical Validation
Constructed mismatch matrices automatically pass the following compliance validations:

| Validation | Description |
|------------|-------------|
| Symmetry | Matrix satisfies `M[i,j] = M[j,i]` |
| Non-negativity | All mismatch counts ≥ 0 |
| Boundedness | All mismatch counts ≤ total loci count |
| Integrality | All mismatch counts are integers |
| Diagonal Rule | Diagonal elements are missing values (NaN) |
| Triangle Inequality | Satisfies `M[i,j] + M[j,k] ≥ M[i,k]` |

### 4. Grading Analysis System
- **Basic Metrics**: Mmin (minimum mismatch), Mmax (maximum mismatch), Mavg (average mismatch), Mzmp (zero-mismatch partner count)
- **Stepwise Classification (Dk)**: Classifies samples into different difference grades based on Mmin values
- **Cumulative Classification (Gk)**: Builds hierarchical system based on cumulative relationship Mmin ≥ k
- **Gap Index (MGI)**: Evaluates the continuity of the grading system
- **Zero-Mismatch Pair Detection**: Identifies Mendelian zero-mismatch sample pairs
- **Duplicate Sample Detection (Mmmd)**: Detects completely duplicate sample groups using Union-Find algorithm

### 5. D0 Network Connectivity Analysis
- Constructs connectivity network for samples with Mmin = 0
- Identifies connected components, isolated samples, and network topology features

### 6. Graphical User Interface (GUI)
- User-friendly interface based on Tkinter
- **Matrix Builder Tab**: Visual input file selection, parameter setting, real-time log viewing
- **Grading Analysis Tab**: Sample filtering, sample comparison, sample search, grading system view, duplicate sample detection, D0 network connectivity analysis

---

## Installation

```bash
pip install numpy pandas matplotlib networkx
```

### System Requirements
- Python ≥ 3.8
- Operating System: Windows / macOS / Linux

---

## Quick Start

### Method 1: Graphical Interface (Recommended for Beginners)

```bash
python MMMv1_en.py
```

Running without parameters automatically launches the GUI:
1. Select **Input CSV File**
2. Select **Output Directory**
3. Click **【Start Build】**

### Method 2: Command Line Mode

#### Build matrix with automatic grading
```bash
python MMMv1_en.py -i input.csv -o output_dir
```

#### Grade existing matrix only
```bash
python MMMv1_en.py --grading -m MMM_Matrix.csv -o grading_output
```

#### Specify SSR marker type
```bash
python MMMv1_en.py -i input.csv -o output_dir -t SSR
```

#### Disable multi-threading (for debugging)
```bash
python MMMv1_en.py -i input.csv -o output_dir --no-threads
```

### Command Line Arguments

| Argument | Description |
|----------|-------------|
| `-i, --input` | Input raw genotype CSV file path |
| `-o, --output` | Output directory path |
| `-m, --matrix` | Existing MMM matrix file path (grading mode) |
| `-t, --type` | Marker type: `SSR` or `SNP` (default: auto-detect) |
| `--no-grading` | Skip automatic grading after matrix construction |
| `--grading` | Execute grading analysis only |
| `--no-threads` | Disable multi-threading |
| `-v, --version` | Show version information |

---

## Input Data Format Examples

### Format 1: Original (Slash-delimited)
```csv
SampleID,LOCUS1,LOCUS2,LOCUS3
Sample001,120/130,140/150,160/170
Sample002,120/120,140/140,160/160
```

### Format 2: Reformatted (Split-column)
```csv
SampleID,LOCUS1_1,LOCUS1_2,LOCUS2_1,LOCUS2_2
Sample001,120,130,140,150
Sample002,120,120,140,140
```

### Format 3: SNP (0/1/2 Coded)
```csv
SampleID,SNP1,SNP2,SNP3
Sample001,0,1,2
Sample002,0,0,1
```

---

## Output Files

After execution, the output directory will contain the following files:

| Filename | Description |
|----------|-------------|
| `MMM_Matrix.csv` | Mismatch matrix M |
| `MMM_Effective_Loci_Matrix.csv` | Effective loci matrix L |
| `MMM_Analytical_Report.csv` | Complete analysis report |
| `MMM_Dk_Classification.csv` | Dk classification results |
| `MMM_basic.csv` | Basic metrics table (Mmin/Mmax/Mavg/Mzmp) |
| `MMM_Grading_Summary.csv` | Grading system summary |
| `MMM_MGI.txt` | Gap index information |
| `MMM_Zero_Mismatch_Pairs.csv` | Zero-mismatch pair detection results |
| `MMM_Duplicate_Samples.csv` | Duplicate sample pairs list |
| `MMM_Duplicate_Groups.csv` | Duplicate sample groups list |
| `MMM_Build_Report.txt` | Build summary report (text format) |

---

## Algorithm Principles

### Mendelian Mismatch Criterion
For two sample genotypes `g1 = (a1, a2)` and `g2 = (b1, b2)`, a Mendelian mismatch is determined when they **share no common alleles**:

```
is_mismatch(g1, g2) = true  if and only if  {a1, a2} ∩ {b1, b2} = ∅
```

### Mismatch Matrix Definition
For a dataset with N samples and L loci, the mismatch matrix M is an N×N symmetric matrix:

```
M[i,j] = Number of loci with Mendelian mismatch between sample i and sample j
```

Diagonal elements `M[i,i]` are set to missing values (NaN).

### Dk Stepwise Classification
The stepwise grade `Dk` of sample `s` is defined as its minimum mismatch count with all other samples:

```
Dk(s) = min{ M[s, j] | j ≠ s }
```

A smaller Dk grade indicates closer genetic relationships within the population.

---

## Project Structure

```
.
├── MMMv1.py          # Main program (Chinese version)
├── MMMv1_en.py       # Main program (English version)
├── README.md         # Chinese documentation
├── README_EN.md      # English documentation
└── ...
```

---

## Version Information

- **Current Version**: v1.0.0
- **Release Date**: 2026-04-19
- **Author**: MMM Assistant

---

## License

This project is open-sourced under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this system in your research, please cite:

> MMM Assistant. (2026). Mendelian Mismatch Matrix Analysis System (MMM v1.0).

---

## Contact

For questions or suggestions, please submit feedback via GitHub Issues.
