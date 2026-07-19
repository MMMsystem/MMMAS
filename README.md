# Mendelian Mismatch Matrix Analysis System (MMMAS)

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

English | [中文](README_ZH.md)

## Introduction

The **Mendelian Mismatch Matrix Analysis System (MMMAS)** is a germplasm bank pre-screening tool based on the Boolean Mendelian exclusion principle. By constructing an N×N mismatch matrix, the system enables duplicate detection, grading evaluation, genetic compatibility analysis, and other functions for germplasm resources.

Key features include:

- Automatic recognition of multiple SSR/SNP genotype data formats
- High-performance calculation of the Mendelian mismatch matrix M, effective loci matrix L, and mismatch rate matrix MMR
- Automatic marker type detection (SSR / SNP)
- Five-fold mathematical compliance validation
- MES (Mendelian Exhaustive Stratification) grading analysis
- Zero-mismatch partner identification and duplicate sample detection
- Complete analytical reports and visualizable outputs
- Both Graphical User Interface (GUI) and command-line operation modes

---

## Core Concepts and Terminology

### Basic Metrics

| Abbreviation | Full Name | Definition |
|---|---|---|
| MMM | Mendelian Mismatch Matrix | N×N matrix where element M[i,j] is the number of Mendelian incompatible loci between sample i and j |
| Mmn | Mendelian Mismatch Number | Number of incompatible loci between two specific samples |
| Mmin | Mendelian Minimum Mismatch Number | Minimum mismatch value of a sample across all pairwise comparisons |
| Mmax | Mendelian Maximum Mismatch Number | Maximum mismatch value of a sample across all pairwise comparisons |
| Mavg | Mendelian Average Mismatch Number | Average mismatch value of a sample across all pairwise comparisons |
| Mzmp | Mendelian Zero-mismatch Partners Number | Number of other samples with zero mismatch (Mmn=0) to a given sample |
| Mmmd | Mendelian Mismatch Mode Duplication | Whether a zero-mismatch pair is a complete duplicate (1=duplicate, 0=non-duplicate) |
| L | Mendelian Effective Loci Matrix | N×N matrix recording the number of effectively comparable loci between each pair of samples |
| MMR | Mendelian Mismatch Rate Matrix | N×N matrix recording the Mendelian mismatch rate between each pair of samples (Mmn/L) |

### Grading System

| Abbreviation | Full Name | Definition |
|---|---|---|
| MES | Mendelian Exhaustive Stratification | Germplasm grading system based on Mmin |
| Dk | Difference Grade k | Set of samples with Mmin = k |
| Ck | Cumulative Grade k | Set of samples with Mmin ≥ k |
| MGI | Mendelian Grade Gap Index | Ratio of empty Dk grades to total grades |

---

## Installation and Runtime Environment

### System Requirements

| Item | Requirement |
|---|---|
| Python | ≥ 3.8 |
| Operating System | Windows / macOS / Linux |
| Memory | ≥ 4GB (8GB or more recommended) |

### Dependencies

```bash
pip install numpy>=1.20.0 pandas>=1.3.0 matplotlib networkx psutil
```

### Installation

Extract the package and run directly from the directory:

```bash
python MMMAS.py
```

---

## Quick Start

### Method 1: Graphical User Interface (Recommended for Beginners)

Running without parameters automatically launches the GUI:

```bash
python MMMAS.py
```

Operation steps:

1. Click "Browse" to select the input CSV file
2. Select the output directory (default: `./MMM_Output/`)
3. Select marker type (optional, auto-detect by default)
4. Set advanced options (optional)
5. Click "Start Analysis"
6. View progress and results

### Method 2: Command Line Mode

#### Basic Usage

```bash
python MMMAS.py -i <input_file> -o <output_dir> [options]
```

#### Common Examples

```bash
# Basic usage
python MMMAS.py -i apple_genotype.csv -o ./apple_results/

# Specify SNP data
python MMMAS.py -i snp_data.csv -m SNP -o ./snp_results/

# Large-scale data (>1000 samples)
python MMMAS.py -i large_dataset.csv -o ./large_results/ -cs 500

# English interface
python MMMAS.py -i data.csv -o ./results/ -l en
```

### Command Line Arguments

| Argument | Short | Description | Default |
|---|---|---|---|
| `--input` | `-i` | Input genotype data file (CSV) | Required |
| `--output` | `-o` | Output directory | `./MMM_Output/` |
| `--marker-type` | `-m` | Marker type (SSR/SNP) | Auto-detect |
| `--no-grading` | `-ng` | Disable automatic grading | False |
| `--no-checkpoint` | `-nc` | Disable checkpoint | False |
| `--no-multithreading` | `-nm` | Disable multi-threading | False |
| `--chunk-size` | `-cs` | Chunk size | 1000 |
| `--language` | `-l` | Interface language (zh/en) | zh |
| `--version` | `-v` | Show version information | - |

---

## Input Data Format

MMMAS supports the following input formats, all automatically detectable:

### Format A: Single-column Allele Split Format (Recommended)

```csv
SampleID,Locus1_1,Locus1_2,Locus2_1,Locus2_2,...
Sample1,100,102,150,154,...
Sample2,100,104,150,150,...
```

### Format B: Single-column Genotype Format

```csv
SampleID,Locus1,Locus2,Locus3,...
Sample1,100/102,150/154,200/204,...
Sample2,100/104,150/150,200/202,...
```

### Format C: SSR Raw Data Format (Double Header)

```csv
SampleID,Locus1,Locus1,Locus2,Locus2,...
,Allele1,Allele2,Allele1,Allele2,...
Sample1,100,102,150,154,...
Sample2,100,104,150,150,...
```

### Format D: SNP Transposed Format

```csv
SampleID,Sample1,Sample2,Sample3,...
Marker1,AA,AT,TT,...
Marker2,CC,CG,GG,...
```

### SNP 0/1/2 Coded Format

```csv
SampleID,SNP1,SNP2,SNP3,...
Sample1,0,1,2
Sample2,1,2,0
```

### Data Requirements

| Item | Requirement |
|---|---|
| Samples | ≥ 2 |
| Loci | ≥ 5 (≥ 11 recommended) |
| Missing data | Allowed, represented by 0 or empty value |
| Encoding | SSR/SNP auto-detection |

---

## Output Files

After execution, the output directory will contain the following files:

### Core Matrix Files

| Filename | Description |
|---|---|
| `MMM_Matrix.csv` | N×N symmetric matrix recording Mendelian incompatible loci between sample i and j; diagonal is NA |
| `MMM_Effective_Loci_Matrix.csv` | N×N matrix recording the number of effectively comparable loci between each pair of samples |
| `MMM_Mmr_Matrix.csv` | N×N matrix recording the Mendelian mismatch rate between each pair of samples (Mmn/L) |

### Grading Files

| Filename | Description |
|---|---|
| `MMM_basic.csv` | Basic metrics per sample: Mmin, Dk Grade, Mmax, Mavg, Mzmp, Mzmp list. |
| `MMM_Grading_Summary.csv` | Sample count, percentage, and member list for each Dk/Ck grade |

### Zero-Mismatch and Duplicate Files

| Filename | Description |
|---|---|
| `MMM_Zero_Mismatch_Pairs.csv` | Zero-mismatch pairs and Mmmd markers |
| `MMM_Duplicate_Samples.csv` | Records with Mmmd=1 only |
| `MMM_Duplicate_Groups.csv` | Duplicate sample groups merged by Union-Find |
| `MMM_only_data.csv` | Genotype data after removing duplicates |

### Frequency Statistics Files

| Filename | Description |
|---|---|
| `MMM_Mnp_Mismatch_Pairs.csv` | List of pairs with mismatch count 0~2 |
| `MMM_frequency.csv` | Frequency distributions of Mmn, Mmin, Mzmp, and Mavg |

### Report Files

| Filename | Description |
|---|---|
| `MMM_Analytical_Report.csv` | Table report with run parameters, validation results, and grading summary |
| `MMM_Build_Report.txt` | Text report with runtime, phase timings, output file list, and descriptions |
| `MMM_MGI.txt` | Mendelian Grade Gap Index, Max Mmin, Total Samples |

---

## Core Module Architecture

```
MMMAS.py
├── MMMConfig              # Configuration management
├── MMMResults             # Results container
├── DataStructureDetector  # Data format auto-detection
├── MarkerTypeDetector     # Marker type auto-detection
├── GenotypeProcessor      # Genotype array construction
├── MismatchCalculator     # Mismatch calculation (standard/vectorized)
├── ParallelMMMCalculator  # Parallel mismatch calculation
├── MMMValidator           # Matrix compliance validation
├── DkClassifier           # MES grading calculation
├── GradingAnalyzer        # Grading analysis (duplicate detection, etc.)
├── ReportGenerator        # Report generation
└── MMMWorkflow            # Unified workflow orchestration
```

Workflow:

```
Input CSV → Format Detection → Genotype Array Construction → Mismatch Matrix Calculation → Matrix Validation → Dk Grading → Zero-Mismatch Detection → Report Generation
```

---

## Algorithm Principles

### Mendelian Incompatibility Criterion

For two diploid genotypes `(a,b)` and `(c,d)`:

- Compatible (M=0) if the two genotypes share at least one allele
- Incompatible (M=1) if the two genotypes share no common alleles

Formal expressions:

- Compatible condition: `{a,b} ∩ {c,d} ≠ ∅`
- Incompatible condition: `{a,b} ∩ {c,d} = ∅`

### Mismatch Count Calculation

For two samples i and j, the mismatch count M[i,j] is the sum of Mendelian incompatibility decisions across all effective loci:

```
M[i,j] = Σ(Incompatible(G[i,l], G[j,l]))  for all effective loci l
```

where `Incompatible()` returns 1 if incompatible and 0 if compatible.

### MES Grading Algorithm

- Calculate Mmin for each sample (row minimum, excluding diagonal)
- `Dk = {samples | Mmin = k}`
- `Ck = {samples | Mmin ≥ k}`
- `MGI = empty Dk grades / (max Mmin + 1)`

### Duplicate Detection Algorithm

1. Find all sample pairs with Mmn=0
2. For each zero-mismatch pair, compare their mismatch profiles with all third-party samples
3. If two samples have exactly the same mismatch profile → Mmmd=1 (duplicate)
4. Use Union-Find to merge duplicate pairs into duplicate groups

---

## Five-Fold Mathematical Validation

Constructed mismatch matrices automatically pass the following compliance validations:

| Validation | Description |
|---|---|
| Symmetry | Matrix satisfies `M[i,j] = M[j,i]` |
| Non-negativity | All mismatch counts ≥ 0 |
| Boundedness | All mismatch counts ≤ total loci count |
| Integrality | All mismatch counts are integers |
| Diagonal Rule | Diagonal elements are missing values (NaN) |

---

## Performance Optimization

### Computation Strategy Selection

| Samples | Strategy | Time Complexity |
|---|---|---|
| ≤ 100 | Standard loop | O(N² × L) |
| ≤ 1,000 | Vectorized computation | O(N² × L) |
| > 1,000 | Parallel block computation | O(N² × L / P) |

### Performance Reference

| Samples | Loci | Estimated Time | Memory |
|---|---|---|---|
| 100 | 15 | < 1 second | ~50MB |
| 500 | 15 | ~5 seconds | ~200MB |
| 1,000 | 15 | ~30 seconds | ~800MB |
| 5,000 | 15 | ~5 minutes | ~4GB |

---

## Frequently Asked Questions

### Q1: How to handle large datasets (>5000 samples)?

Enable multi-threading and adjust the chunk size:

```bash
python MMMAS.py -i large.csv -o ./output/ -cs 500
```

### Q2: What does a high Mzmp value indicate?

A high Mzmp (zero-mismatch partner count) indicates that the sample has multiple genetically identical accessions in the library, possibly:

- A widely cultivated mainstream variety
- A core parent in historical breeding programs
- A synonym (different names for the same variety)

### Q3: What is the difference between Mmmd=0 and Mmmd=1?

- **Mmmd=1**: The two samples have identical mismatch profiles with all other samples, confirming they are duplicates
- **Mmmd=0**: The two samples are zero-mismatch only to each other, but have different mismatch profiles with third-party samples, indicating a very close kinship relationship

### Q4: How to interpret the MGI value?

| MGI Range | Interpretation |
|---|---|
| MGI = 0 | Grades are continuous, no gaps (ideal state) |
| MGI < 0.1 | Grades are basically continuous |
| MGI 0.1~0.3 | Obvious gaps exist |
| MGI > 0.3 | Grades are seriously discontinuous; data quality may need review |

### Q5: What is the biological significance of Dk grading?

| Dk Grade | Biological Significance |
|---|---|
| Dk=0 | Core germplasm, no closer relatives in the library |
| Dk=1 | Has at least one zero-mismatch partner, indicating a close kinship relationship |
| Dk≥2 | Genetic distance gradually increases, kinship relationship becomes more distant |

---

## Version Information

- **Current Version**: v1.0.0
- **Release Date**: 2026-07-19
- **Author**: Qiliang Chen ,Liuxiu Chen

---

## License

This project is open-sourced under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this system in your research, please cite:

[![DOI](https://zenodo.org/badge/1215017107.svg)](https://doi.org/10.5281/zenodo.21441605)

---

## Contact Information

| Item | Information |
|---|---|
| Email | 26451851@qq.com|
| Institution | 湖北省农业科学院果树茶叶研究所|

For questions or suggestions, please contact us through the channels above.
