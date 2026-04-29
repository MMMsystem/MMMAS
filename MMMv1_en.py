#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mendelian Mismatch Matrix Analysis System
Mendelian Mismatch Matrix Analysis System - v1.0

Integrated Features:
  1. Matrix Build: Auto-detect format from raw genotype data and build MMM matrix
  2. Grading Analysis: Based on MMM matrix, calculate Mmin/Mmax/Mavg/Mpc, Dk/Gk classification,
              duplicate sample detection, D0 network connectivity analysis

Author: MMM Assistant
Version: 1.0.0
Date: 2026-04-19
"""

import sys
import os
import re
import time
import argparse
import warnings
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from collections import Counter
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Callable, Union
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

import numpy as np
import pandas as pd

# Network visualization
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

# Version info
__version__ = "4.0.0"
__author__ = "MMM Assistant"

# =============================================================================
# Part 1: Global Configuration and Constants
# =============================================================================

# Missing value definition
MISSING_VALUES = [None, '', 'na', 'nan', 'NA', 'NaN', 'null', 'NULL',
                  '0', './.', '-', '.', '-9', '-1', '000', '999', '9999']

# SNP format special handling: 0 is valid genotype, not missing value
SNP_MISSING_VALUES = [mv for mv in MISSING_VALUES if mv != '0']

# Output filename standardization
OUTPUT_FILES = {
    'mmm_matrix': 'MMM_Matrix.csv',
    'valid_matrix': 'MMM_Effective_Loci_Matrix.csv',
    'report': 'MMM_Analytical_Report.csv',
    'dk_classification': 'MMM_Dk_Classification.csv',
    'summary': 'MMM_Build_Report.txt',
    'basic': 'MMM_basic.csv',
    'grading_summary': 'MMM_Grading_Summary.csv',
    'mgi': 'MMM_MGI.txt',
    'zero_pairs': 'MMM_Zero_Mismatch_Pairs.csv',
    'duplicate_samples': 'MMM_Duplicate_Samples.csv',
    'duplicate_groups': 'MMM_Duplicate_Groups.csv',
}

# Performance config
MAX_WORKERS = min(multiprocessing.cpu_count(), 8)
DEFAULT_CHUNK_SIZE = 500

# GUI constants
GUI_FONT_TITLE = ('Microsoft YaHei', 14, 'bold')
GUI_FONT_LABEL = ('Microsoft YaHei', 9)
GUI_FONT_MONO = ('Consolas', 10)

# Numeric constants
NAN_THRESHOLD = 1e-6
MAX_SMALL_NETWORK = 80
MAX_LARGE_NETWORK = 50

# Column name constants
COLUMN_SAMPLEID = 'SampleID'
COLUMN_MMIN = 'Mmin'
COLUMN_MMAX = 'Mmax'
COLUMN_MAVG = 'Mavg'
COLUMN_MZMP = 'Mzmp'
COLUMN_MMMD = 'Mmmd'


# =============================================================================
# Part 2: Exception Classes
# =============================================================================

class MMMException(Exception):
    """MMM base tool exception"""
    pass

class DataFormatError(MMMException):
    """Data format error"""
    pass

class AnalysisError(MMMException):
    """Analysis execution error"""
    pass

class ExportError(MMMException):
    """Export operation error"""
    pass


# =============================================================================
# Part 3: Data Structure Definitions
# =============================================================================

@dataclass
class MMMConfig:
    """MMM analysis configuration"""
    input_file: str
    output_dir: str
    marker_type: Optional[str] = None
    chunk_size: int = DEFAULT_CHUNK_SIZE
    use_multithreading: bool = True
    enable_checkpoint: bool = True
    min_valid_loci_ratio: float = 0.8
    auto_grading: bool = True


@dataclass
class MMMResults:
    """MMM analysis result container"""
    incomp_mat: np.ndarray = field(default_factory=lambda: np.array([]))
    valid_mat: np.ndarray = field(default_factory=lambda: np.array([]))
    sample_ids: List[str] = field(default_factory=list)
    locus_names: List[str] = field(default_factory=list)
    marker_type: str = ''
    n_samples: int = 0
    n_loci: int = 0
    dk_classification: Dict = field(default_factory=dict)
    validation_report: Dict = field(default_factory=dict)
    quality_report: Optional[Dict] = None
    timestamps: Dict = field(default_factory=dict)
    fingerprint: Dict = field(default_factory=dict)
    grading_result: Dict = field(default_factory=dict)
    zero_mismatch_pairs: Optional[pd.DataFrame] = None


# =============================================================================
# Part 4: Format Detector
# =============================================================================

class DataStructureDetector:
    """Data structure detector - enhanced version"""

    FORMAT_PATTERNS = {
        'slash': re.compile(r'^(\d+)/(\d+)$'),
        'dash': re.compile(r'^(\d+)-(\d+)$'),
        'space': re.compile(r'^(\d+)\s+(\d+)$'),
        'snp_012': re.compile(r'^[012]$'),
    }

    @classmethod
    def detect_format(cls, df: pd.DataFrame) -> str:
        """Detect data format type"""
        if df.shape[1] < 2:
            return 'unknown'

        sample_data = df.iloc[:, 1:].iloc[0].dropna().astype(str).tolist()
        if not sample_data:
            return 'unknown'

        # Check 1: whether slash-delimited
        slash_count = sum(1 for v in sample_data if '/' in str(v))
        if slash_count >= len(sample_data) * 0.5:
            return 'original'

        # Check 2: split-column format
        cols = df.columns.tolist()[1:]
        if len(cols) >= 2:
            paired_count = 0
            for i in range(0, min(len(cols)-1, 20), 2):
                col1, col2 = str(cols[i]), str(cols[i+1])
                if (col1.endswith('_1') and col2.endswith('_2') and col1[:-2] == col2[:-2]) or \
                   (col1.endswith('.1') and col2.endswith('.2') and col1[:-2] == col2[:-2]) or \
                   (col1.endswith('_a') and col2.endswith('_b') and col1[:-2] == col2[:-2]) or \
                   (col2 == col1 + '.1'):
                    paired_count += 1
                elif re.match(r'^.+_[aA]$', col1) and re.match(r'^.+_[bB]$', col2):
                    root1, root2 = col1[:-2], col2[:-2]
                    if root1 == root2:
                        paired_count += 1

            if paired_count >= len(cols) // 4:
                return 'reformatted'

        # Check 3: SNP format (0/1/2)
        all_values = set()
        for col in df.columns[1:]:
            for v in df[col].dropna().astype(str).values[:100]:
                v_clean = v.strip()
                if v_clean.lower() not in ('na', 'nan', '', 'none', 'null', './.', '-', '.'):
                    try:
                        v_int = str(int(float(v_clean)))
                        all_values.add(v_int)
                    except:
                        all_values.add(v_clean)
                        break

        if all_values and all_values.issubset({'0', '1', '2'}):
            return 'snp'

        # Check 4: transposed base SNP format (row=locus, col=sample, value=CC/TC/TT, etc.)
        sample_vals = []
        for col in df.columns[1:min(20, len(df.columns))]:
            sample_vals.extend(df[col].dropna().astype(str).values[:50])

        if sample_vals:
            base_chars = {'A', 'T', 'C', 'G', 'N', 'R', 'Y', 'K', 'M', 'S', 'W', 'H', 'B', 'V', 'D'}
            two_base_count = sum(1 for v in sample_vals if len(str(v).strip()) == 2
                                  and all(ch.upper() in base_chars for ch in str(v).strip()))
            if two_base_count / len(sample_vals) > 0.8:
                return 'transposed_base_snp'

        return 'unknown'

    @classmethod
    def detect_structure(cls, df: pd.DataFrame, format_type: str) -> Dict:
        """Extract data structure info based on format type"""
        sample_ids = df.iloc[:, 0].astype(str).tolist()
        data_cols = df.columns[1:].tolist()

        locus_map = {}
        locus_roots = []

        if format_type == 'reformatted':
            processed = set()
            i = 0
            while i < len(data_cols):
                if i in processed:
                    i += 1
                    continue

                col = data_cols[i]
                col_str = str(col)

                # First detect root/root.1 format
                if i + 1 < len(data_cols):
                    next_col = str(data_cols[i + 1])
                    if next_col == col_str + '.1' or next_col == col_str + '.2':
                        locus_map[col_str] = [col, data_cols[i + 1]]
                        locus_roots.append(col_str)
                        processed.add(i)
                        processed.add(i + 1)
                        i += 2
                        continue

                # Try matching _1/_2, .1/.2, _a/_b patterns
                match = re.match(r'^(.+?)(?:_1|_2|\.1|\.2|_a|_b)$', col_str)
                if match:
                    root = match.group(1)
                    if root not in locus_map:
                        pair_cols = []
                        for suffix in ['_1', '_2']:
                            test_col = root + suffix
                            if test_col in data_cols:
                                pair_cols.append(test_col)
                        if len(pair_cols) == 2:
                            locus_map[root] = pair_cols
                            locus_roots.append(root)
                            processed.add(data_cols.index(pair_cols[0]))
                            processed.add(data_cols.index(pair_cols[1]))
                        else:
                            pair_cols = []
                            for suffix in ['.1', '.2']:
                                test_col = root + suffix
                                if test_col in data_cols:
                                    pair_cols.append(test_col)
                            if len(pair_cols) == 2:
                                locus_map[root] = pair_cols
                                locus_roots.append(root)
                                for pc in pair_cols:
                                    if pc in data_cols:
                                        processed.add(data_cols.index(pc))
                else:
                    # Try explicit a/b suffix
                    match = re.match(r'^(.+?)(?:_a|_b|\.a|\.b)$', col_str)
                    if match and i + 1 < len(data_cols):
                        root = match.group(1)
                        next_col = str(data_cols[i + 1])
                        if next_col == root + '_b' or next_col == root + '_B' or \
                           next_col == root + '.b' or next_col == root + '.B':
                            if root not in locus_map:
                                locus_map[root] = [col, data_cols[i + 1]]
                                locus_roots.append(root)
                                processed.add(i)
                                processed.add(i + 1)
                i += 1

        elif format_type in ('snp', 'original'):
            locus_roots = data_cols
            for col in data_cols:
                locus_map[col] = col

        elif format_type == 'transposed_base_snp':
            # Transposed base SNP format: row=locus, col=sample
            locus_roots = df.iloc[:, 0].astype(str).str.strip().tolist()
            sample_ids = [str(c).strip() for c in data_cols]
            locus_map = {col: col for col in locus_roots}

        return {
            'sample_ids': sample_ids,
            'n_samples': len(sample_ids),
            'locus_map': locus_map,
            'n_loci': len(locus_roots),
            'locus_names': locus_roots,
            'format_type': format_type
        }


# =============================================================================
# Part 5: Marker Type Detector
# =============================================================================

class MarkerTypeDetector:
    """Auto-detect SSR/SNP marker type"""

    @staticmethod
    def detect(df: pd.DataFrame, structure_info: Dict) -> str:
        """Auto-detect marker type"""
        format_type = structure_info.get('format_type', 'unknown')
        if format_type in ('snp', 'transposed_base_snp'):
            return 'SNP'

        all_values = set()
        data_cols = df.columns[1:].tolist()
        sample_cols = data_cols[:min(20, len(data_cols))]

        for col in sample_cols:
            for val in df[col].dropna().values:
                val_str = str(val).strip()
                if '/' in val_str:
                    parts = val_str.split('/')
                    all_values.update(parts)
                else:
                    all_values.add(val_str)

        all_values = {v for v in all_values if v not in MISSING_VALUES}
        is_012 = all_values.issubset({'0', '1', '2'})

        numeric_values = []
        for v in all_values:
            try:
                numeric_values.append(int(float(v)))
            except:
                pass

        if numeric_values:
            max_val = max(numeric_values)
            min_val = min(numeric_values)
            value_range = max_val - min_val
        else:
            max_val = min_val = value_range = 0

        high_loci_count = structure_info.get('n_loci', 0) > 50

        score_snp = 0
        score_ssr = 0

        if is_012:
            score_snp += 3
        if numeric_values and value_range <= 2:
            score_snp += 2
        if high_loci_count:
            score_snp += 1
        if format_type == 'original':
            score_ssr += 1
        if numeric_values and max_val > 100:
            score_ssr += 2

        return 'SNP' if score_snp > score_ssr else 'SSR'


# =============================================================================
# Part 6: Genotype Processor
# =============================================================================

class GenotypeProcessor:
    """Genotype processor - supports SSR and SNP"""

    @staticmethod
    def _standardize_value(val, marker_type: str = 'SSR') -> str:
        """Normalized value is string"""
        if pd.isna(val):
            return ''
        val_str = str(val).strip()
        try:
            fval = float(val_str)
            if marker_type == 'SSR' and fval != int(fval):
                return val_str
            return str(int(fval))
        except (ValueError, TypeError):
            return val_str

    @classmethod
    def build_geno_array(cls, df: pd.DataFrame, structure: Dict, marker_type: str) -> Tuple[np.ndarray, np.ndarray]:
        """Build genotype array"""
        n = structure['n_samples']
        n_loci = structure['n_loci']
        locus_names = structure['locus_names']
        format_type = structure['format_type']

        geno_array = np.full((n, n_loci, 2), None, dtype=object)

        if format_type == 'reformatted':
            for i, loc in enumerate(locus_names):
                cols = structure['locus_map'][loc]
                geno_array[:, i, 0] = df[cols[0]].apply(cls._standardize_value).values
                geno_array[:, i, 1] = df[cols[1]].apply(cls._standardize_value).values

        elif format_type == 'snp':
            for i, loc in enumerate(locus_names):
                col_data = df[loc].apply(lambda x: cls._standardize_value(x, marker_type='SNP')).values
                for j, val in enumerate(col_data):
                    if val == '0':
                        geno_array[j, i] = ['0', '0']
                    elif val == '1':
                        geno_array[j, i] = ['0', '1']
                    elif val == '2':
                        geno_array[j, i] = ['1', '1']
                    else:
                        geno_array[j, i] = [val, val]

        elif format_type == 'transposed_base_snp':
            # Transposed base SNP format: row=locus, col=sample
            # geno_array shape is (n_samples, n_loci, 2)
            # i = sample index (col), j = locus index (row)
            data_cols = df.columns[1:]
            for i, sample_col in enumerate(data_cols):
                col_data = df[sample_col].apply(cls._standardize_value).values
                for j, val in enumerate(col_data):
                    val = str(val).strip()
                    if '/' in val:
                        parts = val.split('/')
                        geno_array[i, j] = [parts[0], parts[1]]
                    elif len(val) == 2:
                        geno_array[i, j] = [val[0], val[1]]
                    elif len(val) == 1:
                        geno_array[i, j] = [val, val]
                    else:
                        geno_array[i, j] = [val, val]

        else:  # original format
            for i, loc in enumerate(locus_names):
                col_data = df[loc].apply(lambda x: cls._standardize_value(x, marker_type=marker_type)).values
                for j, val in enumerate(col_data):
                    if '/' in val:
                        parts = val.split('/')
                        geno_array[j, i] = [parts[0], parts[1]]
                    elif '-' in val:
                        parts = val.split('-')
                        geno_array[j, i] = [parts[0], parts[1]]
                    else:
                        geno_array[j, i] = [val, val]

        # Build missing value mask
        missing_mask = np.zeros((n, n_loci), dtype=bool)
        missing_values = SNP_MISSING_VALUES if marker_type == 'SNP' else MISSING_VALUES

        for mv in missing_values:
            missing_mask |= (geno_array[:, :, 0] == mv) | (geno_array[:, :, 1] == mv)

        missing_mask |= (geno_array[:, :, 0] == '') | (geno_array[:, :, 1] == '')

        return geno_array, missing_mask

    @staticmethod
    def is_mismatch(g1: Tuple, g2: Tuple) -> bool:
        """Determine if two genotypes have Mendelian mismatch"""
        if any(v is None or v == '' for v in [g1[0], g1[1], g2[0], g2[1]]):
            return False

        g1_set = set(g1)
        g2_set = set(g2)
        return len(g1_set & g2_set) == 0



# =============================================================================
# Part 7: Mismatch Matrix Calculator
# =============================================================================

class MismatchCalculator:
    """Mendelian mismatch matrix calculator"""

    def __init__(self, geno_array: np.ndarray, missing_mask: np.ndarray):
        self.geno_array = geno_array
        self.missing_mask = missing_mask
        self.n_samples = geno_array.shape[0]
        self.n_loci = geno_array.shape[1]

    def calculate_standard(self, progress_callback: Optional[Callable] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Standard algorithm (triple loop) - suitable for small datasets"""
        n, n_loci = self.n_samples, self.n_loci
        incomp_mat = np.zeros((n, n), dtype=np.int32)
        valid_mat = np.zeros((n, n), dtype=np.int32)

        total_pairs = n * (n - 1) // 2
        processed = 0

        for i in range(n - 1):
            for j in range(i + 1, n):
                mismatch_count = 0
                valid_count = 0

                for l in range(n_loci):
                    if not self.missing_mask[i, l] and not self.missing_mask[j, l]:
                        valid_count += 1
                        if GenotypeProcessor.is_mismatch(
                            tuple(self.geno_array[i, l]),
                            tuple(self.geno_array[j, l])
                        ):
                            mismatch_count += 1

                incomp_mat[i, j] = incomp_mat[j, i] = mismatch_count
                valid_mat[i, j] = valid_mat[j, i] = valid_count

                processed += 1
                if progress_callback and processed % 1000 == 0:
                    progress_callback(processed / total_pairs * 100)

        return incomp_mat, valid_mat

    def calculate_vectorized(self, progress_callback: Optional[Callable] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Vectorized algorithm - medium-scale datasets"""
        n, n_loci = self.n_samples, self.n_loci
        allele1, allele2 = self._encode_alleles()

        incomp_mat = np.zeros((n, n), dtype=np.int32)
        valid_mat = np.zeros((n, n), dtype=np.int32)

        total_pairs = n * (n - 1) // 2
        processed = 0
        block_size = 128

        for i_start in range(0, n - 1, block_size):
            i_end = min(i_start + block_size, n)

            for j_start in range(i_start + 1, n, block_size):
                j_end = min(j_start + block_size, n)

                a1_i = allele1[i_start:i_end]
                a2_i = allele2[i_start:i_end]
                m_i = self.missing_mask[i_start:i_end]

                a1_j = allele1[j_start:j_end]
                a2_j = allele2[j_start:j_end]
                m_j = self.missing_mask[j_start:j_end]

                ni, nj = i_end - i_start, j_end - j_start

                valid_i = ~m_i
                valid_j = ~m_j
                valid_both = valid_i[:, None, :] & valid_j[None, :, :]
                L_block = np.sum(valid_both, axis=2)

                share1 = (a1_i[:, None, :] == a1_j[None, :, :]) | (a1_i[:, None, :] == a2_j[None, :, :])
                share2 = (a2_i[:, None, :] == a1_j[None, :, :]) | (a2_i[:, None, :] == a2_j[None, :, :])
                has_share = share1 | share2

                M_block = L_block - np.sum(has_share & valid_both, axis=2)

                for ii in range(ni):
                    i_global = i_start + ii
                    for jj in range(nj):
                        j_global = j_start + jj
                        if j_global <= i_global:
                            continue
                        incomp_mat[i_global, j_global] = incomp_mat[j_global, i_global] = int(M_block[ii, jj])
                        valid_mat[i_global, j_global] = valid_mat[j_global, i_global] = int(L_block[ii, jj])
                        processed += 1

                if progress_callback and processed % 1000 == 0:
                    progress_callback(processed / total_pairs * 100)

        return incomp_mat, valid_mat

    @staticmethod
    def _encode_alleles_static(geno_array: np.ndarray, missing_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Encode alleles as integers"""
        n, n_loci = geno_array.shape[0], geno_array.shape[1]
        allele1 = np.full((n, n_loci), -1, dtype=np.int32)
        allele2 = np.full((n, n_loci), -1, dtype=np.int32)

        for l in range(n_loci):
            allele_map = {}
            next_code = 0
            for i in range(n):
                if missing_mask[i, l]:
                    continue
                a1 = geno_array[i, l, 0]
                a2 = geno_array[i, l, 1]
                if a1 not in allele_map:
                    allele_map[a1] = next_code
                    next_code += 1
                allele1[i, l] = allele_map[a1]
                if a2 not in allele_map:
                    allele_map[a2] = next_code
                    next_code += 1
                allele2[i, l] = allele_map[a2]

        return allele1, allele2

    def _encode_alleles(self) -> Tuple[np.ndarray, np.ndarray]:
        return self._encode_alleles_static(self.geno_array, self.missing_mask)


# =============================================================================
# Part 8: High-Performance Parallel Calculator
# =============================================================================

class ParallelMMMCalculator:
    """High-performance parallel MMM matrix calculator"""

    def __init__(self, checkpoint_dir: Optional[str] = None):
        self.checkpoint_dir = checkpoint_dir
        if checkpoint_dir:
            os.makedirs(checkpoint_dir, exist_ok=True)

    def calculate(self, geno_array: np.ndarray, missing_mask: np.ndarray,
                  progress_callback: Optional[Callable] = None,
                  use_multithreading: bool = True,
                  chunk_size: int = DEFAULT_CHUNK_SIZE) -> Tuple[np.ndarray, np.ndarray]:
        """Parallel mismatch matrix calculation"""
        n, n_loci = geno_array.shape[0], geno_array.shape[1]
        allele1, allele2 = self._encode_alleles(geno_array, missing_mask)

        M_matrix = np.full((n, n), np.nan, dtype=np.float32)
        Lij_matrix = np.full((n, n), np.nan, dtype=np.float32)

        tasks = []
        for i_start in range(0, n, chunk_size):
            i_end = min(i_start + chunk_size, n)
            for j_start in range(i_start, n, chunk_size):
                j_end = min(j_start + chunk_size, n)
                tasks.append((allele1, allele2, missing_mask, i_start, i_end, j_start, j_end))

        total_tasks = len(tasks)
        completed = 0

        if use_multithreading and total_tasks > 1:
            with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
                futures = {executor.submit(self._compute_block, task): task for task in tasks}
                for future in as_completed(futures):
                    (i_start, j_start), M_block, L_block = future.result()
                    i_end = i_start + M_block.shape[0]
                    j_end = j_start + M_block.shape[1]

                    M_matrix[i_start:i_end, j_start:j_end] = M_block
                    Lij_matrix[i_start:i_end, j_start:j_end] = L_block

                    if i_start != j_start:
                        M_matrix[j_start:j_end, i_start:i_end] = M_block.T
                        Lij_matrix[j_start:j_end, i_start:i_end] = L_block.T

                    completed += 1
                    if progress_callback:
                        progress_callback(completed / total_tasks * 100)
        else:
            for task in tasks:
                (i_start, j_start), M_block, L_block = self._compute_block(task)
                i_end = i_start + M_block.shape[0]
                j_end = j_start + L_block.shape[1]

                M_matrix[i_start:i_end, j_start:j_end] = M_block
                Lij_matrix[i_start:i_end, j_start:j_end] = L_block

                if i_start != j_start:
                    M_matrix[j_start:j_end, i_start:i_end] = M_block.T
                    Lij_matrix[j_start:j_end, i_start:i_end] = L_block.T

                completed += 1
                if progress_callback:
                    progress_callback(completed / total_tasks * 100)

        np.fill_diagonal(M_matrix, np.nan)
        np.fill_diagonal(Lij_matrix, np.nan)

        return M_matrix, Lij_matrix

    def _encode_alleles(self, geno_array: np.ndarray, missing_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        return MismatchCalculator._encode_alleles_static(geno_array, missing_mask)

    def _compute_block(self, args) -> Tuple[Tuple[int, int], np.ndarray, np.ndarray]:
        """Compute matrix block"""
        allele1, allele2, missing_mask, i_start, i_end, j_start, j_end = args

        a1_i = allele1[i_start:i_end]
        a2_i = allele2[i_start:i_end]
        a1_j = allele1[j_start:j_end]
        a2_j = allele2[j_start:j_end]
        m_i = missing_mask[i_start:i_end]
        m_j = missing_mask[j_start:j_end]

        ni, nj = i_end - i_start, j_end - j_start
        n_loci = allele1.shape[1]

        M_block = np.zeros((ni, nj), dtype=np.float32)
        L_block = np.zeros((ni, nj), dtype=np.float32)

        for l in range(n_loci):
            valid_i = ~m_i[:, l]
            valid_j = ~m_j[:, l]
            valid_both = valid_i[:, None] & valid_j[None, :]

            share1 = (a1_i[:, l:l+1] == a1_j[None, :, l]) | (a1_i[:, l:l+1] == a2_j[None, :, l])
            share2 = (a2_i[:, l:l+1] == a1_j[None, :, l]) | (a2_i[:, l:l+1] == a2_j[None, :, l])
            has_share = share1 | share2

            mismatch = valid_both & ~has_share

            M_block += mismatch.astype(np.float32)
            L_block += valid_both.astype(np.float32)

        return (i_start, j_start), M_block, L_block


# =============================================================================
# Part 9: Matrix Validator
# =============================================================================

class MMMValidator:
    """MMM matrix compliance validator"""

    def __init__(self, incomp_mat: np.ndarray, valid_mat: np.ndarray,
                 n_loci: int, sample_ids: List[str]):
        self.incomp_mat = incomp_mat
        self.valid_mat = valid_mat
        self.n_loci = n_loci
        self.sample_ids = sample_ids
        self.n_samples = len(sample_ids)

    def validate(self) -> Dict:
        """Execute complete validation process"""
        report = {
            'validation_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'tests': {},
            'overall_status': 'PENDING'
        }

        tests = [
            ('symmetry', 'Symmetry Validation', self._test_symmetry),
            ('nonnegative', 'Non-negativity Validation', self._test_nonnegative),
            ('bounded', 'Boundedness Validation', self._test_bounded),
            ('integer', 'Integrality Validation', self._test_integer),
            ('diagonal', 'Diagonal Rule Validation', self._test_diagonal),
            ('triangle_inequality', 'Triangle Inequality Validation', self._test_triangle),
        ]

        passed = 0
        for test_id, test_name, test_func in tests:
            try:
                result = test_func()
                report['tests'][test_id] = {
                    'name': test_name,
                    'status': 'PASS' if result else 'FAIL',
                    'requirement': '-',
                }
                if result:
                    passed += 1
            except Exception as e:
                report['tests'][test_id] = {
                    'name': test_name,
                    'status': f'Error: {str(e)}',
                    'requirement': '-',
                }

        report['passed_count'] = passed
        report['test_count'] = len(tests)
        report['compliant'] = passed == len(tests)
        report['overall_status'] = 'ALL PASSED' if report['compliant'] else 'PARTIALLY FAILED'

        return report

    def _test_symmetry(self) -> bool:
        m = np.nan_to_num(self.incomp_mat, nan=0)
        return np.allclose(m, m.T)

    def _test_nonnegative(self) -> bool:
        """Non-negativity validation: all non-missing values >= 0"""
        valid = self.incomp_mat[~np.isnan(self.incomp_mat)]
        return len(valid) > 0 and np.all(valid >= 0)

    def _test_bounded(self) -> bool:
        """Boundedness validation: all non-missing values <= n_loci"""
        valid = self.incomp_mat[~np.isnan(self.incomp_mat)]
        return len(valid) > 0 and np.all(valid <= self.n_loci)

    def _test_integer(self) -> bool:
        """Integrality validation: all non-missing values are integers"""
        valid = self.incomp_mat[~np.isnan(self.incomp_mat)]
        return len(valid) > 0 and np.all(np.equal(np.mod(valid, 1), 0))

    def _test_diagonal(self) -> bool:
        """Diagonal rule validation: diagonal elements are NaN"""
        return np.all(np.isnan(np.diag(self.incomp_mat)))

    def _test_triangle(self) -> bool:
        """Triangle inequality validation"""
        for i in range(self.n_samples):
            for j in range(i + 1, self.n_samples):
                for k in range(j + 1, self.n_samples):
                    m_ij = self.incomp_mat[i, j]
                    m_jk = self.incomp_mat[j, k]
                    m_ik = self.incomp_mat[i, k]
                    if not any(np.isnan([m_ij, m_jk, m_ik])):
                        if m_ij + m_jk < m_ik or m_ij + m_ik < m_jk or m_jk + m_ik < m_ij:
                            return False
        return True


# =============================================================================
# Part 10: Dk Classifier
# =============================================================================

class DkClassifier:
    """MES classification calculator"""

    def __init__(self, incomp_mat: np.ndarray, sample_ids: List[str]):
        self.incomp_mat = incomp_mat
        self.sample_ids = sample_ids
        self.n_samples = len(sample_ids)

    def classify(self) -> Dict:
        """Execute Dk classification"""
        result = {
            'mmin': {},
            'mmax': {},
            'mmean': {},
            'mzmp': {},
            'dk': {},
            'dk_groups': {},
            'gk_groups': {},
            'gap_index': 0.0,
            'mzmp_partners': {}
        }

        for i, sid in enumerate(self.sample_ids):
            row = self.incomp_mat[i, :]
            valid_mask = ~np.isnan(row)

            if not np.any(valid_mask):
                continue

            valid_values = row[valid_mask]

            mmin = int(np.min(valid_values))
            result['mmin'][sid] = mmin

            mmax = int(np.max(valid_values))
            result['mmax'][sid] = mmax

            mmean = float(np.mean(valid_values))
            result['mmean'][sid] = mmean

            mzmp = int(np.sum(valid_values == 0))
            result['mzmp'][sid] = mzmp

            partners = [self.sample_ids[j] for j in range(self.n_samples)
                        if i != j and not np.isnan(row[j]) and row[j] == 0]
            result['mzmp_partners'][sid] = partners

            result['dk'][sid] = mmin

        max_k = max(result['mmin'].values()) if result['mmin'] else 0
        for k in range(max_k + 1):
            result['dk_groups'][k] = [
                sid for sid, v in result['mmin'].items() if v == k
            ]

        for k in range(max_k + 1):
            result['gk_groups'][k] = [
                sid for sid, v in result['mmin'].items() if v >= k
            ]

        empty_grades = sum(1 for k in range(max_k + 1)
                          if k not in result['dk_groups'] or len(result['dk_groups'][k]) == 0)
        result['gap_index'] = empty_grades / (max_k + 1) if max_k > 0 else 0

        return result




# =============================================================================
# Part 11: Grading Analyzer (integrate MatrixGradingAnalyzer)
# =============================================================================

class GradingAnalyzer:
    """
    Matrix grading analysis core class
    Integrates core features of MMM_Matrix_Grading.py
    """

    @staticmethod
    def validate_matrix(mmm_matrix_df: pd.DataFrame) -> None:
        """Validate input matrix format and content"""
        if mmm_matrix_df is None or mmm_matrix_df.empty:
            raise DataFormatError("Input matrix is empty")

        if mmm_matrix_df.shape[0] < 2:
            raise DataFormatError("Too few samples in matrix (need at least 2)")

        if mmm_matrix_df.shape[1] < 2:
            raise DataFormatError("Too few columns in matrix (need ID column + at least 1 sample column)")

        try:
            matrix_values = mmm_matrix_df.iloc[:, 1:].values
            matrix_values.astype(float)
        except (ValueError, TypeError) as e:
            raise DataFormatError(f"Matrix data columns contain non-numeric data: {str(e)}")

    @staticmethod
    def analyze(mmm_matrix_df: pd.DataFrame) -> Dict:
        """
        Perform comprehensive analysis based on mismatch matrix

        Returns a dict with the following keys:
            - basic_df: Basic metrics table (Mmin/Mmax/Mavg/Mpc)
            - grading_summary_df: Grading system summary table
            - zero_pairs_df: Zero-mismatch pair info table
            - duplicate_pairs_df: Duplicate sample pair table
            - duplicate_groups: Duplicate sample group list
            - duplicate_sample_list: Sample IDs involved in duplicates
            - G_dict: Cumulative grade dictionary
            - D_dict: Stepwise grade dictionary
            - MGI: Gap index
            - max_mmin: Max Mmin value
            - sample_list: Sample ID list
            - mmin_list: Mmin list
            - M_matrix: Mismatch matrix
        """
        try:
            GradingAnalyzer.validate_matrix(mmm_matrix_df)

            sample_id_col = mmm_matrix_df.columns[0]
            sample_list = mmm_matrix_df[sample_id_col].tolist()
            N = len(sample_list)

            M_matrix = mmm_matrix_df.iloc[:, 1:].values.astype(float)

            # ===================== Basic Metrics Calculation =====================
            M_diag_masked = np.ma.masked_array(M_matrix, mask=np.eye(N, dtype=bool))
            mmin_list = np.ma.getdata(np.ma.min(M_diag_masked, axis=1)).astype(int)
            mmax_list = np.ma.getdata(np.ma.max(M_diag_masked, axis=1)).astype(int)

            mavg_list = []
            mpc_list = []
            for i in range(N):
                row = np.concatenate([M_matrix[i, :i], M_matrix[i, i+1:]])
                mavg_list.append(round(np.nanmean(row), 2))
                mpc_list.append(int(np.nansum(row == 0)))

            mavg_list = np.array(mavg_list)
            mpc_list = np.array(mpc_list)

            partners_list = []
            for i in range(N):
                partners = [sample_list[j] for j in range(N)
                            if i != j and not np.isnan(M_matrix[i, j]) and M_matrix[i, j] == 0]
                partners_list.append(','.join(partners) if partners else '')

            basic_df = pd.DataFrame({
                sample_id_col: sample_list,
                COLUMN_MMIN: mmin_list,
                'DkLevel': mmin_list,
                COLUMN_MMAX: mmax_list,
                COLUMN_MAVG: mavg_list,
                COLUMN_MZMP: mpc_list,
                'Mzmp list.': partners_list,
            })

            # ===================== Zero-Mismatch Pairs & Duplicate Sample Detection =====================
            zero_pairs = []
            duplicate_pairs = []

            zero_mask = (M_matrix == 0)
            np.fill_diagonal(zero_mask, False)
            i_idx, j_idx = np.where(np.triu(zero_mask, k=1))

            for i, j in zip(i_idx, j_idx):
                mask = np.ones(N, dtype=bool)
                mask[i] = False
                mask[j] = False

                row_i = M_matrix[i, mask]
                row_j = M_matrix[j, mask]
                is_duplicate = np.allclose(row_i, row_j, atol=NAN_THRESHOLD)

                mmmd_val = 1 if is_duplicate else 0
                zero_pairs.append({
                    'SampleID_1': sample_list[i],
                    'SampleID_2': sample_list[j],
                    COLUMN_MMMD: mmmd_val
                })

                if is_duplicate:
                    duplicate_pairs.append({
                        'SampleID_1': sample_list[i],
                        'SampleID_2': sample_list[j],
                        COLUMN_MMMD: 1
                    })

            zero_pairs_df = pd.DataFrame(zero_pairs) if zero_pairs else pd.DataFrame(
                columns=['SampleID_1', 'SampleID_2', COLUMN_MMMD])
            duplicate_pairs_df = pd.DataFrame(duplicate_pairs) if duplicate_pairs else pd.DataFrame(
                columns=['SampleID_1', 'SampleID_2', COLUMN_MMMD])

            # ===================== Duplicate Sample Groups (Union-Find) =====================
            duplicate_groups, duplicate_sample_list = GradingAnalyzer._find_duplicate_groups(
                duplicate_pairs, sample_list
            )

            # ===================== Grading System =====================
            max_mmin = int(mmin_list.max())
            n = max_mmin

            G_dict = {}
            for k in range(0, n + 1):
                members = [sample_list[i] for i in range(N) if mmin_list[i] >= k]
                G_dict[k] = members

            D_dict = {}
            for k in range(0, n + 1):
                members = [sample_list[i] for i in range(N) if mmin_list[i] == k]
                D_dict[k] = members

            empty_D_count = sum(1 for k in range(0, n + 1) if len(D_dict[k]) == 0)
            mgi = empty_D_count / (n + 1) if (n + 1) > 0 else 0.0

            # Grading summary table
            grading_summary = []
            for k in range(0, n + 1):
                grading_summary.append({
                    'Grade_k': k,
                    'G_k_Count': len(G_dict[k]),
                    'G_k_Percent': round(len(G_dict[k]) / N * 100, 2),
                    'G_k_Members': ';'.join(G_dict[k]),
                    'D_k_Count': len(D_dict[k]),
                    'D_k_Percent': round(len(D_dict[k]) / N * 100, 2),
                    'D_k_Members': ';'.join(D_dict[k]),
                    'D_k_Empty': 'Yes' if len(D_dict[k]) == 0 else 'No'
                })
            grading_summary_df = pd.DataFrame(grading_summary)

            return {
                'basic_df': basic_df,
                'grading_summary_df': grading_summary_df,
                'zero_pairs_df': zero_pairs_df,
                'duplicate_pairs_df': duplicate_pairs_df,
                'duplicate_groups': duplicate_groups,
                'duplicate_sample_list': duplicate_sample_list,
                'G_dict': G_dict,
                'D_dict': D_dict,
                'MGI': mgi,
                'max_mmin': n,
                'sample_list': sample_list,
                'mmin_list': mmin_list.tolist(),
                'M_matrix': M_matrix,
            }

        except DataFormatError:
            raise
        except Exception as e:
            raise AnalysisError(f"Analysis execution failed: {str(e)}")

    @staticmethod
    def _find_duplicate_groups(duplicate_pairs: list, sample_list: list) -> tuple:
        """Use union-find to merge duplicate sample pairs into duplicate sample groups"""
        if not duplicate_pairs:
            return [], []

        sample_to_idx = {s: i for i, s in enumerate(sample_list)}
        n = len(sample_list)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        for pair in duplicate_pairs:
            i = sample_to_idx.get(pair['SampleID_1'])
            j = sample_to_idx.get(pair['SampleID_2'])
            if i is not None and j is not None:
                union(i, j)

        groups = {}
        for i in range(n):
            root = find(i)
            groups.setdefault(root, []).append(sample_list[i])

        duplicate_groups = [sorted(g) for g in groups.values() if len(g) > 1]
        duplicate_sample_list = sorted({s for g in duplicate_groups for s in g})

        return duplicate_groups, duplicate_sample_list


# =============================================================================
# Part 12: Report Generator
# =============================================================================

class ReportGenerator:
    """Unified report generator"""

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def save_matrices(self, results: MMMResults) -> None:
        """Save mismatch matrix and effective loci matrix"""
        incomp_df = pd.DataFrame(
            results.incomp_mat,
            index=results.sample_ids,
            columns=results.sample_ids
        )
        incomp_df.insert(0, 'SampleID', results.sample_ids)
        incomp_path = os.path.join(self.output_dir, OUTPUT_FILES['mmm_matrix'])
        incomp_df.to_csv(incomp_path, index=False, encoding='utf-8-sig', na_rep='NA')

        valid_df = pd.DataFrame(
            results.valid_mat,
            index=results.sample_ids,
            columns=results.sample_ids
        )
        valid_df.insert(0, 'SampleID', results.sample_ids)
        valid_path = os.path.join(self.output_dir, OUTPUT_FILES['valid_matrix'])
        valid_df.to_csv(valid_path, index=False, encoding='utf-8-sig', na_rep='NA')

        print(f"  ✓ Mismatch matrix M: {incomp_path}")
        print(f"  ✓ Effective loci matrix L: {valid_path}")

    def save_dk_classification(self, dk_results: Dict, sample_ids: List[str]) -> None:
        """Save Dk classification results"""
        dk_data = []
        for sid in sample_ids:
            if sid in dk_results['mmin']:
                partners = dk_results.get('mzmp_partners', {}).get(sid, [])
                partners_str = ','.join(partners) if partners else ''
                dk_data.append({
                    'SampleID': sid,
                    'MinMismatch': dk_results['mmin'][sid],
                    'MaxMismatch(Dk)': dk_results['mmax'][sid],
                    'AvgMismatch': round(dk_results['mmean'][sid], 2),
                    'ZeroMismatchPartners': dk_results['mzmp'][sid],
                    'ZeroMismatchPartnerList': partners_str
                })

        dk_df = pd.DataFrame(dk_data)
        dk_path = os.path.join(self.output_dir, OUTPUT_FILES['dk_classification'])
        dk_df.to_csv(dk_path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Dk classification table: {dk_path}")

    def save_zero_mismatch_pairs(self, zero_pairs_df: pd.DataFrame) -> None:
        """Save zero-mismatch pairs and Mmmd results"""
        path = os.path.join(self.output_dir, OUTPUT_FILES['zero_pairs'])
        zero_pairs_df.to_csv(path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Zero-mismatch pair table: {path}")
        if not zero_pairs_df.empty:
            dup_count = int((zero_pairs_df['Mmmd'] == 1).sum())
            print(f"      Total {len(zero_pairs_df)} zero-mismatch pairs, including {dup_count} complete duplicates (Mmmd=1)")

    def save_grading_results(self, grading: Dict) -> None:
        """Save grading analysis results"""
        # MMM_basic.csv
        basic_path = os.path.join(self.output_dir, OUTPUT_FILES['basic'])
        grading['basic_df'].to_csv(basic_path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Basic metrics table: {basic_path}")

        # MMM_Grading_Summary.csv
        summary_path = os.path.join(self.output_dir, OUTPUT_FILES['grading_summary'])
        grading['grading_summary_df'].to_csv(summary_path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Grading summary table: {summary_path}")

        # MMM_MGI.txt
        mgi_path = os.path.join(self.output_dir, OUTPUT_FILES['mgi'])
        with open(mgi_path, 'w', encoding='utf-8') as f:
            f.write(f"MGI (gap index) = {grading['MGI']:.4f}\n")
            f.write(f"Max Mmin = {grading['max_mmin']}\n")
            f.write(f"Total samples = {len(grading['sample_list'])}\n")
        print(f"  ✓ MGI index: {mgi_path}")

        # MMM_Zero_Mismatch_Pairs.csv
        zp_path = os.path.join(self.output_dir, OUTPUT_FILES['zero_pairs'])
        grading['zero_pairs_df'].to_csv(zp_path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Zero-mismatch pair table: {zp_path}")

        # MMM_Duplicate_Samples.csv
        if not grading['duplicate_pairs_df'].empty:
            dup_path = os.path.join(self.output_dir, OUTPUT_FILES['duplicate_samples'])
            grading['duplicate_pairs_df'].to_csv(dup_path, index=False, encoding='utf-8-sig')
            print(f"  ✓ Duplicate sample pairs: {dup_path}")

        # MMM_Duplicate_Groups.csv
        if grading['duplicate_groups']:
            groups_data = []
            for idx, group in enumerate(grading['duplicate_groups'], 1):
                groups_data.append({
                    'Group_ID': idx,
                    'Group_Size': len(group),
                    'Members': ';'.join(group)
                })
            dg_path = os.path.join(self.output_dir, OUTPUT_FILES['duplicate_groups'])
            pd.DataFrame(groups_data).to_csv(dg_path, index=False, encoding='utf-8-sig')
            print(f"  ✓ Duplicate sample groups: {dg_path}")

    def generate_analytical_report(self, results: MMMResults,
                                    validation_report: Dict,
                                    total_time: float) -> pd.DataFrame:
        """Generate complete analysis report"""
        rows = []
        n, n_loci = results.n_samples, results.n_loci

        rows.extend([
            ['Basic Info', 'Total Samples (N)', str(n), '≥1', 'Pass' if n >= 1 else 'Insufficient samples'],
            ['Basic Info', 'Total Loci (L)', str(n_loci), '≥5', 'Pass' if n_loci >= 5 else 'Insufficient loci'],
            ['Basic Info', 'Marker Type', results.marker_type, 'SSR/SNP', 'Detected'],
            ['Basic Info', 'Total Build Time (s)', f"{total_time:.2f}", '-', 'Complete'],
        ])

        rows.extend([
            ['Algorithm Info', 'Core Algorithm', 'Vectorized parallel computation', '-', 'MMM standard compliant'],
            ['Algorithm Info', 'Time Complexity', 'O(N²)', '-', 'Optimized'],
        ])

        for test_name, test_info in validation_report['tests'].items():
            rows.append([
                'Mathematical Validation',
                test_info['name'],
                test_info['status'],
                test_info['requirement'],
                test_info['status']
            ])

        rows.append([
            'Mathematical Validation',
            'Final Validation Conclusion',
            validation_report['overall_status'],
            'ALL PASSED',
            validation_report['overall_status']
        ])

        if results.grading_result:
            grading = results.grading_result
            rows.append(['Grading Analysis', 'Gap Index (MGI)', f"{grading['MGI']:.4f}", '-', '-'])
            rows.append(['Grading Analysis', 'Max Mmin', str(grading['max_mmin']), '-', '-'])
            rows.append(['Grading Analysis', 'Zero-Mismatch Pairs', str(len(grading['zero_pairs_df'])), '-', '-'])
            rows.append(['Grading Analysis', 'Duplicate Sample Pairs', str(len(grading['duplicate_pairs_df'])), '-', '-'])
            rows.append(['Grading Analysis', 'Duplicate Sample Groups', str(len(grading['duplicate_groups'])), '-', '-'])

        report_df = pd.DataFrame(rows, columns=['Category', 'Metric', 'Value', 'Standard', 'Conclusion'])
        report_path = os.path.join(self.output_dir, OUTPUT_FILES['report'])
        report_df.to_csv(report_path, index=False, encoding='utf-8-sig')
        print(f"  ✓ Analysis report: {report_path}")

        return report_df

    def generate_summary_txt(self, results: MMMResults, validation_report: Dict,
                             total_time: float) -> str:
        """Generate summary text report"""
        lines = [
            "=" * 80,
            "          Mendelian Mismatch Matrix System (MMM) Analysis Report",
            "=" * 80,
            "",
            f"Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Marker type: {results.marker_type}",
            f"Sample count: {results.n_samples}",
            f"Locus count: {results.n_loci}",
            f"Total time: {total_time:.2f}s",
            "",
            "-" * 80,
            "[FIVE-FOLD MATHEMATICAL VALIDATION]",
            "-" * 80,
        ]

        for test_name, test_info in validation_report['tests'].items():
            icon = "✓" if test_info['status'] == 'PASS' else "✗"
            lines.append(f"  {icon} {test_info['name']}: {test_info['status']}")

        lines.extend([
            "",
            f"  Final conclusion: {validation_report['overall_status']}",
            "",
            "=" * 80,
            "Output files:",
            "=" * 80,
        ])

        for key, filename in OUTPUT_FILES.items():
            lines.append(f"  • {filename}")

        lines.append("=" * 80)

        summary = "\n".join(lines)
        summary_path = os.path.join(self.output_dir, OUTPUT_FILES['summary'])
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(summary)
        print(f"  ✓ Summary report: {summary_path}")

        return summary




# =============================================================================
# Part 13: Unified Workflow
# =============================================================================

class MMMWorkflow:
    """MMM unified workflow: matrix build + grading analysis"""

    def __init__(self, config: MMMConfig):
        self.config = config
        self.df = None
        self.structure = None
        self.marker_type = None
        self.geno_array = None
        self.missing_mask = None

    def run(self) -> MMMResults:
        """Execute complete workflow"""
        start_time = time.time()
        results = MMMResults()
        timings = {}

        # ========== Phase 1: Data Loading ==========
        print("\n[Phase 1] Data loading and format detection...")
        phase1_start = time.time()
        self._load_data()
        phase1_time = time.time() - phase1_start
        timings['Phase 1: Data Loading and Format Detection'] = phase1_time
        print(f"  ⏱ Time: {phase1_time:.2f}s")

        # ========== Phase 2: Build Genotype Array ==========
        print("\n[Phase 2] Build genotype database...")
        phase2_start = time.time()
        self.geno_array, self.missing_mask = GenotypeProcessor.build_geno_array(
            self.df, self.structure, self.marker_type
        )
        phase2_time = time.time() - phase2_start
        timings['Phase 2: Build Genotype Database'] = phase2_time
        print(f"  ✓ Genotype matrix: {self.structure['n_samples']}×{self.structure['n_loci']}×2")
        print(f"  ⏱ Time: {phase2_time:.2f}s")

        # ========== Phase 3: Calculate Mismatch Matrix ==========
        print("\n[Phase 3] Calculate Mendelian mismatch matrix...")
        phase3_start = time.time()
        incomp_mat, valid_mat = self._calculate_mismatch()
        phase3_time = time.time() - phase3_start
        timings['Phase 3: Calculate Mendelian Mismatch Matrix'] = phase3_time
        print(f"  ⏱ Time: {phase3_time:.2f}s")

        # ========== Phase 4: Mathematical Validation ==========
        print("\n[Phase 4] Execute six-fold compliance validation...")
        phase4_start = time.time()
        validator = MMMValidator(incomp_mat, valid_mat, self.structure['n_loci'], self.structure['sample_ids'])
        validation_report = validator.validate()
        phase4_time = time.time() - phase4_start
        timings['Phase 4: Six-fold Compliance Validation'] = phase4_time

        for test_name, test_info in validation_report['tests'].items():
            icon = "✓" if test_info['status'] == 'PASS' else "✗"
            print(f"  {icon} {test_info['name']}: {test_info['status']}")
        print(f"\n  Final conclusion: {validation_report['overall_status']}")
        print(f"  ⏱ Time: {phase4_time:.2f}s")

        # ========== Phase 5: Dk Classification ==========
        print("\n[Phase 5] Calculate Dk classification...")
        phase5_start = time.time()
        classifier = DkClassifier(incomp_mat, self.structure['sample_ids'])
        dk_results = classifier.classify()
        phase5_time = time.time() - phase5_start
        timings['Phase 5: Calculate Dk Classification'] = phase5_time
        print(f"  ✓ Mmin range: 0-{max(dk_results['mmin'].values()) if dk_results['mmin'] else 0}")
        print(f"  ✓ Gap index: {dk_results['gap_index']:.3f}")
        print(f"  ⏱ Time: {phase5_time:.2f}s")

        # ========== Phase 5b: Zero-Mismatch Pairs ==========
        print("\n[Phase 5b] Detect zero-mismatch pairs and Mmmd...")
        phase5b_start = time.time()
        zero_pairs_df = self._find_zero_mismatch_pairs(incomp_mat)
        phase5b_time = time.time() - phase5b_start
        timings['Phase 5b: Detect Zero-Mismatch Pairs and Mmmd'] = phase5b_time
        print(f"  ✓ Zero-mismatch pairs: {len(zero_pairs_df)} pairs")
        if not zero_pairs_df.empty:
            dup_count = int((zero_pairs_df['Mmmd'] == 1).sum())
            print(f"  ✓ Complete duplicates (Mmmd=1): {dup_count} pairs")
        print(f"  ⏱ Time: {phase5b_time:.2f}s")

        # ========== Phase 6: Grading Analysis (Optional) ==========
        grading_result = {}
        phase6_time = 0
        if self.config.auto_grading:
            print("\n[Phase 6] Execute grading analysis...")
            phase6_start = time.time()
            incomp_df = pd.DataFrame(
                incomp_mat,
                index=self.structure['sample_ids'],
                columns=self.structure['sample_ids']
            )
            incomp_df.insert(0, 'SampleID', self.structure['sample_ids'])
            grading_result = GradingAnalyzer.analyze(incomp_df)
            phase6_time = time.time() - phase6_start
            timings['Phase 6: Grading Analysis'] = phase6_time
            print(f"  ✓ Grading system: 0-{grading_result['max_mmin']} levels")
            print(f"  ✓ MGI: {grading_result['MGI']:.4f}")
            print(f"  ✓ Duplicate sample groups: {len(grading_result['duplicate_groups'])} groups")
            print(f"  ⏱ Time: {phase6_time:.2f}s")

        end_time = time.time()
        total_time = end_time - start_time

        # Assemble result object for report generation
        results = MMMResults(
            incomp_mat=incomp_mat,
            valid_mat=valid_mat,
            sample_ids=self.structure['sample_ids'],
            locus_names=self.structure['locus_names'],
            marker_type=self.marker_type,
            n_samples=self.structure['n_samples'],
            n_loci=self.structure['n_loci'],
            dk_classification=dk_results,
            validation_report=validation_report,
            grading_result=grading_result,
            zero_mismatch_pairs=zero_pairs_df,
            timestamps={
                'analysis_start': datetime.fromtimestamp(start_time).isoformat(),
                'analysis_end': datetime.fromtimestamp(end_time).isoformat(),
                'total_seconds': round(total_time, 3),
                'phase_timings': {k: round(v, 3) for k, v in timings.items()},
            }
        )

        # ========== Phase 7: Generate Reports ==========
        print("\n[Phase 7] Generate analysis report...")
        phase7_start = time.time()
        generator = ReportGenerator(self.config.output_dir)
        generator.save_matrices(results)
        generator.save_zero_mismatch_pairs(zero_pairs_df)
        if grading_result:
            generator.save_grading_results(grading_result)
        generator.generate_analytical_report(results, validation_report, total_time)
        generator.generate_summary_txt(results, validation_report, total_time)
        phase7_time = time.time() - phase7_start
        timings['Phase 7: Generate Analysis Reports'] = phase7_time
        print(f"  ⏱ Time: {phase7_time:.2f}s")

        print("\n" + "=" * 80)
        print("Phase Timing Summary:")
        for phase_name, seconds in timings.items():
            print(f"  {phase_name}: {seconds:.2f}s")
        print(f"{'-'*40}")
        print(f"Total time: {total_time:.2f}s")
        print(f"Output directory: {self.config.output_dir}")
        print("=" * 80)

        return results

    def run_grading_only(self, matrix_file: str) -> Dict:
        """Execute grading analysis only (from existing matrix file)"""
        print("=" * 80)
        print("       MMM Matrix Grading Analysis Mode")
        print("=" * 80)
        print(f"\nLoading matrix file: {matrix_file}")

        mmm_df = pd.read_csv(matrix_file, encoding='utf-8-sig')
        N = mmm_df.shape[0]
        print(f"  Sample count: {N}")
        print(f"  Matrix dimensions: {N} x {N}")

        print("\n[Step 1] Calculate basic metrics (Mmin/Mmax/Mavg/Mpc)...")
        grading_result = GradingAnalyzer.analyze(mmm_df)
        print(f"  ✓ Max Mmin: {grading_result['max_mmin']}")
        print(f"  ✓ MGI (gap index): {grading_result['MGI']:.4f}")

        print(f"\n[Step 2] Detect zero-mismatch pairs and duplicate samples...")
        zero_count = len(grading_result['zero_pairs_df'])
        dup_count = len(grading_result['duplicate_pairs_df'])
        group_count = len(grading_result['duplicate_groups'])
        print(f"  ✓ Zero-mismatch pairs: {zero_count} pairs")
        print(f"  ✓ Complete duplicates (Mmmd=1): {dup_count} pairs")
        print(f"  ✓ Duplicate sample groups: {group_count} groups")

        print(f"\n[Step 3] Build grading system (Dk/Gk)...")
        for k in range(0, grading_result['max_mmin'] + 1):
            d_count = len(grading_result['D_dict'][k])
            g_count = len(grading_result['G_dict'][k])
            print(f"  D{k}: {d_count:3d} items  |  G{k}: {g_count:3d} items")

        print(f"\n[Step 4] Export results...")
        generator = ReportGenerator(self.config.output_dir)
        generator.save_grading_results(grading_result)

        # Also save MGI
        mgi_path = os.path.join(self.config.output_dir, OUTPUT_FILES['mgi'])
        with open(mgi_path, 'w', encoding='utf-8') as f:
            f.write(f"MGI (gap index) = {grading_result['MGI']:.4f}\n")
            f.write(f"Max Mmin = {grading_result['max_mmin']}\n")
            f.write(f"Total samples = {len(grading_result['sample_list'])}\n")
        print(f"  ✓ MGI index: {mgi_path}")

        print("\n" + "=" * 80)
        print("Grading analysis complete!")
        print(f"  Output directory: {self.config.output_dir}")
        print("=" * 80)

        return grading_result

    def _load_data(self):
        """Load data"""
        encodings = ['utf-8', 'utf-8-sig', 'gbk', 'gb2312', 'latin-1']
        read_success = False
        last_error = None

        for encoding in encodings:
            try:
                self.df = pd.read_csv(self.config.input_file, low_memory=False, dtype=str, encoding=encoding)
                print(f"  ✓ File encoding: {encoding}")
                read_success = True
                break
            except Exception as e:
                last_error = e
                continue

        if not read_success:
            raise ValueError(f"Cannot read file, please ensure it is a valid CSV format. Last error: {last_error}")

        format_type = DataStructureDetector.detect_format(self.df)
        if format_type == 'unknown':
            raise ValueError("Cannot auto-detect data format, please check file format")

        print(f"  ✓ Data format: {format_type}")

        self.structure = DataStructureDetector.detect_structure(self.df, format_type)
        print(f"  ✓ Samples: {self.structure['n_samples']}")
        print(f"  ✓ Loci: {self.structure['n_loci']}")

        if self.config.marker_type:
            self.marker_type = self.config.marker_type
            print(f"  ✓ Marker type (user-specified): {self.marker_type}")
        else:
            self.marker_type = MarkerTypeDetector.detect(self.df, self.structure)
            print(f"  ✓ Marker type (auto-detected): {self.marker_type}")

    def _calculate_mismatch(self) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate mismatch matrix"""
        n = self.structure['n_samples']

        if n <= 100:
            calculator = MismatchCalculator(self.geno_array, self.missing_mask)
            incomp_mat, valid_mat = calculator.calculate_standard()
        elif n <= 1000:
            calculator = MismatchCalculator(self.geno_array, self.missing_mask)
            incomp_mat, valid_mat = calculator.calculate_vectorized()
        else:
            calc = ParallelMMMCalculator(
                checkpoint_dir=self.config.output_dir if self.config.enable_checkpoint else None
            )
            incomp_mat, valid_mat = calc.calculate(
                self.geno_array, self.missing_mask,
                use_multithreading=self.config.use_multithreading,
                chunk_size=self.config.chunk_size
            )

        incomp_mat = incomp_mat.astype(float)
        valid_mat = valid_mat.astype(float)
        np.fill_diagonal(incomp_mat, np.nan)
        np.fill_diagonal(valid_mat, np.nan)

        total_pairs = n * (n - 1) // 2
        print(f"  ✓ Completed {total_pairs} pair calculations")

        return incomp_mat, valid_mat

    def _find_zero_mismatch_pairs(self, incomp_mat: np.ndarray) -> pd.DataFrame:
        """Detect zero-mismatch pairs and Mmmd"""
        sample_ids = self.structure['sample_ids']
        n = len(sample_ids)
        zero_pairs = []

        zero_mask = (incomp_mat == 0)
        np.fill_diagonal(zero_mask, False)
        i_idx, j_idx = np.where(np.triu(zero_mask, k=1))

        for i, j in zip(i_idx, j_idx):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            mask[j] = False

            row_i = incomp_mat[i, mask]
            row_j = incomp_mat[j, mask]
            is_duplicate = np.allclose(row_i, row_j, atol=NAN_THRESHOLD)

            mmmd_val = 1 if is_duplicate else 0
            zero_pairs.append({
                'SampleID_1': sample_ids[i],
                'SampleID_2': sample_ids[j],
                'Mmmd': mmmd_val
            })

        df = pd.DataFrame(zero_pairs) if zero_pairs else pd.DataFrame(
            columns=['SampleID_1', 'SampleID_2', 'Mmmd'])
        return df




# =============================================================================
# Part 14: Graphical User Interface
# =============================================================================

class ThreadSafeTextRedirector:
    """Thread-safe text redirector"""

    def __init__(self, queue_obj: queue.Queue):
        self.queue = queue_obj

    def write(self, s: str):
        if s:
            self.queue.put(s)

    def flush(self):
        pass

    def isatty(self):
        return False


class MMMApp:
    """MMM integrated graphical user interface"""

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Mendelian Mismatch Matrix Analysis System - v1.0")
        self.root.geometry("980x760")
        self.root.minsize(800, 600)

        self.font_ui = ("Microsoft YaHei", 10)
        self.font_title = ("Microsoft YaHei", 14, "bold")
        self.font_mono = ("Consolas", 10)

        # Output queue
        self.output_queue: queue.Queue = queue.Queue()
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.redirector = ThreadSafeTextRedirector(self.output_queue)

        self._build_ui()
        self._poll_queue()

    def _build_ui(self):
        """Build interface"""
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(
            title_frame,
            text="Mendelian Mismatch Matrix Analysis System",
            font=self.font_title,
            anchor='center'
        ).pack(fill=tk.X)
        ttk.Label(
            title_frame,
            text="Mendelian Mismatch Matrix Analysis System - v1.0",
            font=("Microsoft YaHei", 10),
            anchor='center'
        ).pack(fill=tk.X)

        # Notebook
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1: Matrix Build
        self.builder_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.builder_tab, text="  Matrix Build  ")
        self._build_builder_tab()

        # Tab 2: Grading Analysis
        self.grading_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(self.grading_tab, text="  Grading Analysis  ")
        self._build_grading_tab()

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(
            main_frame, textvariable=self.status_var,
            font=("Microsoft YaHei", 9, "italic"), foreground="gray"
        ).pack(fill=tk.X, pady=(5, 0))

        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    # ==================== Matrix Build Tab ====================

    def _build_builder_tab(self):
        """Build matrix build tab"""
        frame = self.builder_tab
        frame.columnconfigure(1, weight=1)

        # Input file
        ttk.Label(frame, text="Input CSV File:", font=self.font_ui).grid(row=0, column=0, sticky=tk.W, pady=5)
        self.builder_input_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.builder_input_var, font=self.font_mono).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(frame, text="Browse...", command=self._browse_builder_input).grid(row=0, column=2, pady=5)

        # Output directory
        ttk.Label(frame, text="Output Directory:", font=self.font_ui).grid(row=1, column=0, sticky=tk.W, pady=5)
        self.builder_output_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.builder_output_var, font=self.font_mono).grid(
            row=1, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(frame, text="Browse...", command=self._browse_builder_output).grid(row=1, column=2, pady=5)

        # Marker type
        ttk.Label(frame, text="Marker Type:", font=self.font_ui).grid(row=2, column=0, sticky=tk.W, pady=5)
        self.marker_var = tk.StringVar(value="Auto Detect")
        ttk.Combobox(frame, textvariable=self.marker_var, values=["Auto Detect", "SSR", "SNP"],
                     state="readonly", width=15).grid(row=2, column=1, sticky=tk.W, padx=5, pady=5)

        # Advanced options
        options_frame = ttk.LabelFrame(frame, text="Advanced Options", padding="10")
        options_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        self.threads_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Enable Multithreading", variable=self.threads_var).grid(row=0, column=0, sticky=tk.W, padx=5)

        self.auto_grading_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Auto-run grading analysis after build", variable=self.auto_grading_var).grid(row=0, column=1, sticky=tk.W, padx=20)

        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=4, column=0, columnspan=3, pady=10)

        self.run_builder_btn = ttk.Button(btn_frame, text="▶ Start Build", command=self._start_builder, width=18)
        self.run_builder_btn.pack(side=tk.LEFT, padx=5)

        # Instructions
        help_frame = ttk.LabelFrame(frame, text="Instructions", padding="5")
        help_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        help_text = (
            "1. Select raw genotype CSV file → 2. Select output directory → 3. Click [Start Build]\n"
            "Supported formats: original(slash-delimited) / reformatted(split-column) / SNP(0/1/2 encoded)\n"
            "After build, auto-fill MMM_Matrix.csv path into [Grading Analysis] tab"
        )
        ttk.Label(help_frame, text=help_text, font=("Microsoft YaHei", 8), foreground="gray").pack(anchor=tk.W)

        # Log
        log_frame = ttk.LabelFrame(frame, text="Build Log", padding="5")
        log_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        frame.rowconfigure(6, weight=1)

        self.builder_log = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, font=self.font_mono,
                                                      state=tk.DISABLED, height=18)
        self.builder_log.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    def _browse_builder_input(self):
        filename = filedialog.askopenfilename(title="Select Input CSV File",
                                               filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if filename:
            self.builder_input_var.set(filename)
            if not self.builder_output_var.get():
                default_out = os.path.join(os.path.dirname(filename), "MMM_Output")
                self.builder_output_var.set(default_out)

    def _browse_builder_output(self):
        dirname = filedialog.askdirectory(title="Select Output Directory")
        if dirname:
            self.builder_output_var.set(dirname)

    def _start_builder(self):
        input_file = self.builder_input_var.get().strip()
        output_dir = self.builder_output_var.get().strip()

        if not input_file or not os.path.isfile(input_file):
            messagebox.showerror("Error", "Please select a valid input CSV file")
            return
        if not output_dir:
            messagebox.showerror("Error", "Please select an output directory")
            return

        self.run_builder_btn.configure(state=tk.DISABLED)
        self.status_var.set("Building matrix, please wait...")
        self.builder_log.configure(state=tk.NORMAL)
        self.builder_log.delete(1.0, tk.END)
        self.builder_log.configure(state=tk.DISABLED)

        sys.stdout = self.redirector
        sys.stderr = self.redirector

        marker_type = self.marker_var.get()
        if marker_type == "Auto Detect":
            marker_type = None

        config = MMMConfig(
            input_file=input_file,
            output_dir=output_dir,
            marker_type=marker_type,
            use_multithreading=self.threads_var.get(),
            auto_grading=self.auto_grading_var.get()
        )

        thread = threading.Thread(target=self._run_builder_thread, args=(config,), daemon=True)
        thread.start()

    def _run_builder_thread(self, config: MMMConfig):
        try:
            workflow = MMMWorkflow(config)
            results = workflow.run()
            self.root.after(0, self._on_builder_success, results)
        except Exception as e:
            self.root.after(0, lambda err=e: self._on_builder_error(err))
        finally:
            sys.stdout = self.original_stdout
            sys.stderr = self.original_stderr
            self.root.after(0, self._reset_builder_ui)

    def _on_builder_success(self, results: MMMResults):
        msg = "Matrix build completed successfully!"
        if results.grading_result:
            msg += f"\n\nGrading analysis result:\n  MGI: {results.grading_result['MGI']:.4f}\n  Duplicate sample groups: {len(results.grading_result['duplicate_groups'])} groups"

        # Auto-fill generated matrix path into grading analysis input
        output_dir = self.builder_output_var.get().strip()
        if output_dir:
            matrix_path = os.path.join(output_dir, OUTPUT_FILES['mmm_matrix'])
            if os.path.isfile(matrix_path):
                self.grading_input_var.set(matrix_path)
                # Also set grading analysis output directory to Grading subdir under matrix output dir
                grading_out = os.path.join(output_dir, "Grading_Output")
                self.grading_output_var.set(grading_out)
                msg += f"\n\nAuto-filled grading analysis matrix file:\n  {matrix_path}"

        messagebox.showinfo("Analysis Complete", msg)

    def _on_builder_error(self, error: Exception):
        messagebox.showerror("Analysis Failed", f"❌ Error during analysis:\n\n{str(error)}")

    def _reset_builder_ui(self):
        self.run_builder_btn.configure(state=tk.NORMAL)
        self.status_var.set("Ready")

    # ==================== Grading Analysis Tab ====================

    def _build_grading_tab(self):
        """Build grading analysis tab"""
        frame = self.grading_tab
        frame.columnconfigure(1, weight=1)

        # Input matrix
        ttk.Label(frame, text="MMM Matrix File:", font=self.font_ui).grid(row=0, column=0, sticky=tk.W, pady=5)
        self.grading_input_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.grading_input_var, font=self.font_mono).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(frame, text="Browse...", command=self._browse_grading_input).grid(row=0, column=2, pady=5)

        # Output directory
        ttk.Label(frame, text="Output Directory:", font=self.font_ui).grid(row=1, column=0, sticky=tk.W, pady=5)
        self.grading_output_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.grading_output_var, font=self.font_mono).grid(
            row=1, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(frame, text="Browse...", command=self._browse_grading_output).grid(row=1, column=2, pady=5)

        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=2, column=0, columnspan=3, pady=10)

        buttons = [
            ("Start Analysis", self._start_grading),
            ("Sample View & Filter", self._view_sample_filter),
            ("Sample Compare", self._view_sample_comparison),
            ("Sample Search", self._view_sample_search),
            ("Grading System", self._view_grading_system),
            ("Duplicate Samples", self._view_duplicate_samples),
            ("D0 Network", self._view_d0_network),
            ("Clear Log", self._clear_grading_log),
        ]

        for label, command in buttons:
            ttk.Button(btn_frame, text=label, command=command, width=12).pack(side=tk.LEFT, padx=3)

        # Instructions
        help_frame = ttk.LabelFrame(frame, text="Instructions", padding="5")
        help_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        help_text = (
            "1. Select MMM_Matrix.csv file → 2. Select output directory → 3. Click [Start Analysis]\n"
            "After analysis, use: Sample View & Filter / Sample Compare / Sample Search / Grading System / Duplicate Samples / D0 Network"
        )
        ttk.Label(help_frame, text=help_text, font=("Microsoft YaHei", 8), foreground="gray").pack(anchor=tk.W)

        # Log
        log_frame = ttk.LabelFrame(frame, text="Analysis Log", padding="5")
        log_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        frame.rowconfigure(4, weight=1)

        self.grading_log = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, font=self.font_mono,
                                                      state=tk.DISABLED, height=20)
        self.grading_log.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        self.grading_result = None

    def _browse_grading_input(self):
        filename = filedialog.askopenfilename(title="Select MMM_Matrix.csv",
                                               filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if filename:
            self.grading_input_var.set(filename)
            if not self.grading_output_var.get():
                self.grading_output_var.set(os.path.join(os.path.dirname(filename), "Grading_Output"))

    def _browse_grading_output(self):
        dirname = filedialog.askdirectory(title="Select Output Directory")
        if dirname:
            self.grading_output_var.set(dirname)

    def _start_grading(self):
        input_file = self.grading_input_var.get().strip()
        output_dir = self.grading_output_var.get().strip()

        if not input_file or not os.path.isfile(input_file):
            messagebox.showerror("Error", "Please select a valid MMM matrix file")
            return
        if not output_dir:
            messagebox.showerror("Error", "Please select an output directory")
            return

        self.status_var.set("Performing grading analysis, please wait...")
        self.grading_log.configure(state=tk.NORMAL)
        self.grading_log.delete(1.0, tk.END)
        self.grading_log.configure(state=tk.DISABLED)

        sys.stdout = self.redirector
        sys.stderr = self.redirector

        config = MMMConfig(
            input_file=input_file,
            output_dir=output_dir,
            auto_grading=False
        )

        thread = threading.Thread(target=self._run_grading_thread, args=(config,), daemon=True)
        thread.start()

    def _run_grading_thread(self, config: MMMConfig):
        try:
            workflow = MMMWorkflow(config)
            result = workflow.run_grading_only(config.input_file)
            self.root.after(0, self._on_grading_success, result)
        except Exception as e:
            self.root.after(0, lambda err=e: self._on_grading_error(err))
        finally:
            sys.stdout = self.original_stdout
            sys.stderr = self.original_stderr
            self.root.after(0, lambda: self.status_var.set("Ready"))

    def _on_grading_success(self, result: Dict):
        self.grading_result = result
        messagebox.showinfo("Analysis Complete",
                            f"✅ Grading analysis complete!\n\nMGI: {result['MGI']:.4f}\n"
                            f"Duplicate sample groups: {len(result['duplicate_groups'])} groups\n"
                            f"Output directory: {self.grading_output_var.get()}")

    def _on_grading_error(self, error: Exception):
        messagebox.showerror("Analysis Failed", f"❌ Grading analysis failed:\n\n{str(error)}")

    def _fill_sample_detail(self, sid: str, text_widget: scrolledtext.ScrolledText):
        """Fill detailed sample info into the text widget"""
        try:
            sample_list = self.grading_result['sample_list']
            idx = sample_list.index(sid)
            row = self.grading_result['basic_df'].iloc[idx]
            mmin = int(row[COLUMN_MMIN])

            G_ks = [k for k, members in self.grading_result['G_dict'].items() if sid in members]
            D_k = mmin

            dup_info = ""
            for grp_idx, group in enumerate(self.grading_result['duplicate_groups'], 1):
                if sid in group:
                    others = [g for g in group if g != sid]
                    dup_info = f"\n{'='*50}\n⚠ Duplicate sample reminder: This sample forms a duplicate group with {', '.join(others)} (Group {grp_idx})\n{'='*50}"
                    break

            lines = [
                f"Sample ID: {sid}",
                "=" * 50,
                f"Mmin (minimum mismatch): {row[COLUMN_MMIN]}",
                f"Mmax (maximum mismatch): {row[COLUMN_MMAX]}",
                f"Mavg (average mismatch): {row[COLUMN_MAVG]}",
                f"Mzmp (Mendelian zero-mismatch partners): {row[COLUMN_MZMP]}",
                "=" * 50,
                f"Stepwise D_k: D{D_k}",
                f"Cumulative G_k includes: {sorted(G_ks)}",
                "=" * 50,
                f"Same stepwise D{D_k} members ({len(self.grading_result['D_dict'][D_k])} items):",
            ]
            for m in self.grading_result['D_dict'][D_k]:
                lines.append(f"  - {m}")
            lines.append("=" * 50)
            if dup_info:
                lines.append(dup_info)

            text_widget.delete(1.0, tk.END)
            text_widget.insert(tk.END, "\n".join(lines))
        except Exception as e:
            messagebox.showerror("Error", f"Error updating details: {str(e)}")

    def _view_sample_filter(self):
        """Sample filter - filter by Mmin/Mmax/Mavg range"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        win = tk.Toplevel(self.root)
        win.title("Sample Filter")
        win.geometry("700x550")
        win.transient(self.root)
        win.grab_set()

        basic_df = self.grading_result['basic_df']
        max_mmin = self.grading_result['max_mmin']
        max_mmax = int(basic_df[COLUMN_MMAX].max()) if len(basic_df) > 0 else 0
        max_mavg = float(basic_df[COLUMN_MAVG].max()) if len(basic_df) > 0 else 0.0

        # Filter conditions
        filter_frame = ttk.LabelFrame(win, text="Filter Conditions", padding=5)
        filter_frame.pack(fill=tk.X, padx=10, pady=5)

        # Mmin
        ttk.Label(filter_frame, text="Mmin:").grid(row=0, column=0, sticky=tk.W, padx=2)
        mmin_min_var = tk.IntVar(value=0)
        mmin_max_var = tk.IntVar(value=max_mmin)
        ttk.Spinbox(filter_frame, from_=0, to=max_mmin, textvariable=mmin_min_var, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(filter_frame, text="-").grid(row=0, column=2)
        ttk.Spinbox(filter_frame, from_=0, to=max_mmin, textvariable=mmin_max_var, width=5).grid(row=0, column=3, padx=2)

        # Mmax
        ttk.Label(filter_frame, text="Mmax:").grid(row=1, column=0, sticky=tk.W, padx=2, pady=5)
        mmax_min_var = tk.IntVar(value=0)
        mmax_max_var = tk.IntVar(value=max_mmax)
        ttk.Spinbox(filter_frame, from_=0, to=max_mmax, textvariable=mmax_min_var, width=5).grid(row=1, column=1, padx=2)
        ttk.Label(filter_frame, text="-").grid(row=1, column=2)
        ttk.Spinbox(filter_frame, from_=0, to=max_mmax, textvariable=mmax_max_var, width=5).grid(row=1, column=3, padx=2)

        # Mavg
        ttk.Label(filter_frame, text="Mavg:").grid(row=2, column=0, sticky=tk.W, padx=2)
        mavg_min_var = tk.DoubleVar(value=0.0)
        mavg_max_var = tk.DoubleVar(value=max_mavg)
        ttk.Entry(filter_frame, textvariable=mavg_min_var, width=7).grid(row=2, column=1, padx=2)
        ttk.Label(filter_frame, text="-").grid(row=2, column=2)
        ttk.Entry(filter_frame, textvariable=mavg_max_var, width=7).grid(row=2, column=3, padx=2)

        # Mzmp
        max_mzmp = int(basic_df[COLUMN_MZMP].max()) if len(basic_df) > 0 else 0
        ttk.Label(filter_frame, text="Mzmp:").grid(row=3, column=0, sticky=tk.W, padx=2, pady=5)
        mzmp_min_var = tk.IntVar(value=0)
        mzmp_max_var = tk.IntVar(value=max_mzmp)
        ttk.Spinbox(filter_frame, from_=0, to=max_mzmp, textvariable=mzmp_min_var, width=5).grid(row=3, column=1, padx=2)
        ttk.Label(filter_frame, text="-").grid(row=3, column=2)
        ttk.Spinbox(filter_frame, from_=0, to=max_mzmp, textvariable=mzmp_max_var, width=5).grid(row=3, column=3, padx=2)

        # Result list
        cols = ("SampleID", "Mmin", "Mmax", "Mavg", "Mzmp")
        tree = ttk.Treeview(win, columns=cols, show="headings", height=12)
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=110, anchor='center')
        tree.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Detail text
        detail_text = scrolledtext.ScrolledText(win, wrap=tk.WORD, font=self.font_mono, height=8)
        detail_text.pack(fill=tk.X, padx=10, pady=5)

        def apply_filter():
            for item in tree.get_children():
                tree.delete(item)

            mmin_lo = mmin_min_var.get()
            mmin_hi = mmin_max_var.get()
            mmax_lo = mmax_min_var.get()
            mmax_hi = mmax_max_var.get()
            mavg_lo = mavg_min_var.get()
            mavg_hi = mavg_max_var.get()
            mzmp_lo = mzmp_min_var.get()
            mzmp_hi = mzmp_max_var.get()

            for idx, sid in enumerate(self.grading_result['sample_list']):
                row = basic_df.iloc[idx]
                mmin = float(row[COLUMN_MMIN])
                mmax = float(row[COLUMN_MMAX])
                mavg = float(row[COLUMN_MAVG])
                mzmp = int(row[COLUMN_MZMP])

                if (mmin_lo <= mmin <= mmin_hi and
                    mmax_lo <= mmax <= mmax_hi and
                    mavg_lo <= mavg <= mavg_hi and
                    mzmp_lo <= mzmp <= mzmp_hi):
                    tree.insert("", tk.END, values=(sid, int(mmin), int(mmax), round(mavg, 2), mzmp))

        def on_select(event):
            sel = tree.selection()
            if not sel:
                return
            sid = tree.item(sel[0], 'values')[0]
            self._fill_sample_detail(sid, detail_text)

        tree.bind('<<TreeviewSelect>>', on_select)

        ttk.Button(filter_frame, text="Filter", command=apply_filter).grid(row=0, column=4, rowspan=4, padx=10)
        apply_filter()

        # Export button (top-right)
        def export_filter_results():
            rows = []
            for item in tree.get_children():
                vals = tree.item(item, 'values')
                rows.append({
                    'SampleID': vals[0],
                    'Mmin': vals[1],
                    'Mmax': vals[2],
                    'Mavg': vals[3],
                    'Mzmp': vals[4],
                })
            if not rows:
                messagebox.showwarning("Warning", "No filter results to export")
                return
            df = pd.DataFrame(rows)
            path = filedialog.asksaveasfilename(
                title="Save Filter Results",
                defaultextension=".csv",
                filetypes=[("CSV Files", "*.csv")]
            )
            if path:
                df.to_csv(path, index=False, encoding='utf-8-sig')
                messagebox.showinfo("Export Complete", f"Saved: {path}")

        ttk.Button(filter_frame, text="Export CSV", command=export_filter_results).grid(row=0, column=5, rowspan=4, padx=10)
        ttk.Button(win, text="Close", command=win.destroy).pack(side=tk.BOTTOM, pady=5)

    def _view_sample_comparison(self):
        """Sample comparator - compare mismatch metrics between two samples"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        win = tk.Toplevel(self.root)
        win.title("Sample Comparator")
        win.geometry("500x550")
        win.transient(self.root)
        win.grab_set()

        sample_list = self.grading_result['sample_list']
        basic_df = self.grading_result['basic_df']
        M_matrix = self.grading_result['M_matrix']

        ttk.Label(win, text="Sample A:").pack(pady=2)
        var_a = tk.StringVar()
        combo_a = ttk.Combobox(win, textvariable=var_a, values=sample_list, state='readonly', width=40)
        combo_a.pack(pady=2)

        ttk.Label(win, text="Sample B:").pack(pady=2)
        var_b = tk.StringVar()
        combo_b = ttk.Combobox(win, textvariable=var_b, values=sample_list, state='readonly', width=40)
        combo_b.pack(pady=2)

        result_text = scrolledtext.ScrolledText(win, wrap=tk.WORD, font=self.font_mono, height=22)
        result_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        def do_compare():
            a = var_a.get()
            b = var_b.get()
            if not a or not b:
                messagebox.showwarning("Warning", "Please select two samples")
                return
            if a == b:
                messagebox.showwarning("Warning", "Please select two different samples")
                return

            try:
                idx_a = sample_list.index(a)
                idx_b = sample_list.index(b)
                row_a = basic_df.iloc[idx_a]
                row_b = basic_df.iloc[idx_b]
                m_ab = int(M_matrix[idx_a, idx_b])

                lines = [
                    f"{'='*50}",
                    f"Sample comparison: {a}  vs  {b}",
                    f"{'='*50}",
                    "",
                    f"[{a}]",
                    f"  Mmin: {int(row_a[COLUMN_MMIN])}",
                    f"  Mmax: {int(row_a[COLUMN_MMAX])}",
                    f"  Mavg: {row_a[COLUMN_MAVG]:.2f}",
                    f"  Mzmp: {int(row_a[COLUMN_MZMP])}",
                    "",
                    f"[{b}]",
                    f"  Mmin: {int(row_b[COLUMN_MMIN])}",
                    f"  Mmax: {int(row_b[COLUMN_MMAX])}",
                    f"  Mavg: {row_b[COLUMN_MAVG]:.2f}",
                    f"  Mzmp: {int(row_b[COLUMN_MZMP])}",
                    "",
                    f"{'='*50}",
                    f"A vs B mismatch count M: {m_ab}",
                    f"{'='*50}",
                ]

                # Dk info
                D_a = int(row_a[COLUMN_MMIN])
                D_b = int(row_b[COLUMN_MMIN])
                lines.append(f"")
                lines.append(f"A stepwise D{D_a}: {len(self.grading_result['D_dict'].get(D_a, []))} samples total")
                lines.append(f"B stepwise D{D_b}: {len(self.grading_result['D_dict'].get(D_b, []))} samples total")

                # Duplicate info
                dup_a = ""
                for grp_idx, group in enumerate(self.grading_result['duplicate_groups'], 1):
                    if a in group:
                        others = [g for g in group if g != a]
                        dup_a = f"  ⚠ A belongs to duplicate group {grp_idx}: {', '.join(others)}"
                        break
                dup_b = ""
                for grp_idx, group in enumerate(self.grading_result['duplicate_groups'], 1):
                    if b in group:
                        others = [g for g in group if g != b]
                        dup_b = f"  ⚠ B belongs to duplicate group {grp_idx}: {', '.join(others)}"
                        break
                if dup_a or dup_b:
                    lines.append(f"")
                    lines.append(f"[Duplicate Sample Info]")
                    if dup_a:
                        lines.append(dup_a)
                    if dup_b:
                        lines.append(dup_b)

                result_text.delete(1.0, tk.END)
                result_text.insert(tk.END, "\n".join(lines))
            except Exception as e:
                messagebox.showerror("Error", f"Error during comparison: {str(e)}")

        ttk.Button(win, text="Compare", command=do_compare).pack(pady=2)
        ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)

    def _view_sample_search(self):
        """Sample search locator - search samples by name"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        win = tk.Toplevel(self.root)
        win.title("Sample Search")
        win.geometry("600x550")
        win.transient(self.root)
        win.grab_set()

        ttk.Label(win, text="Search keyword (partial match supported):").pack(pady=2)
        search_var = tk.StringVar()
        ttk.Entry(win, textvariable=search_var, width=45).pack(pady=2)

        cols = ("SampleID", "Mmin", "Mmax", "Mavg", "Mzmp")
        tree = ttk.Treeview(win, columns=cols, show="headings", height=14)
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=100, anchor='center')
        tree.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        detail_text = scrolledtext.ScrolledText(win, wrap=tk.WORD, font=self.font_mono, height=8)
        detail_text.pack(fill=tk.X, padx=10, pady=5)

        def do_search():
            keyword = search_var.get().strip().lower()
            for item in tree.get_children():
                tree.delete(item)

            for idx, sid in enumerate(self.grading_result['sample_list']):
                if not keyword or keyword in sid.lower():
                    row = self.grading_result['basic_df'].iloc[idx]
                    tree.insert("", tk.END, values=(
                        sid,
                        int(row[COLUMN_MMIN]),
                        int(row[COLUMN_MMAX]),
                        round(float(row[COLUMN_MAVG]), 2),
                        int(row[COLUMN_MZMP])
                    ))

        def on_select(event):
            sel = tree.selection()
            if not sel:
                return
            sid = tree.item(sel[0], 'values')[0]
            self._fill_sample_detail(sid, detail_text)

        ttk.Button(win, text="Search", command=do_search).pack(pady=2)
        tree.bind('<<TreeviewSelect>>', on_select)
        ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)
        do_search()  # Default: show all

    def _view_grading_system(self):
        """Grading system and duplicate samples (merged view)"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        win = tk.Toplevel(self.root)
        win.title("Grading System and Duplicate Samples")
        win.geometry("700x750")
        win.transient(self.root)
        win.grab_set()

        notebook = ttk.Notebook(win)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # ========== Tab 1: Grading System ==========
        tab1 = ttk.Frame(notebook)
        notebook.add(tab1, text="Grading System")

        ctrl_frame = ttk.Frame(tab1)
        ctrl_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(ctrl_frame, text="Select Grade_k:").pack(side=tk.LEFT, padx=2)
        k_var = tk.IntVar(value=0)
        max_k = self.grading_result['max_mmin']
        combo = ttk.Combobox(ctrl_frame, textvariable=k_var, values=list(range(0, max_k + 1)),
                             state='readonly', width=8)
        combo.pack(side=tk.LEFT, padx=2)
        combo.current(0)

        text1 = scrolledtext.ScrolledText(tab1, wrap=tk.WORD, font=self.font_mono, width=80, height=36)
        text1.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        def refresh():
            try:
                k = k_var.get()
                N = len(self.grading_result['sample_list'])
                G_members = self.grading_result['G_dict'].get(k, [])
                D_members = self.grading_result['D_dict'].get(k, [])
                G_next = self.grading_result['G_dict'].get(k + 1, [])

                lines = [
                    f"Grade_k = {k}",
                    "=" * 60,
                    f"Cumulative G_{k}: {len(G_members)} items  ({len(G_members)/N*100:.2f}%)",
                ]
                for m in G_members:
                    lines.append(f"  - {m}")
                lines.append("-" * 60)
                lines.append(f"Stepwise D_{k}: {len(D_members)} items  ({len(D_members)/N*100:.2f}%)")
                for m in D_members:
                    lines.append(f"  - {m}")
                lines.append("-" * 60)
                if k < max_k:
                    lines.append(f"Next cumulative G_{k+1}: {len(G_next)} items")
                    lines.append(f"Validation: D_{k} = G_{k} \\ G_{{{k+1}}}  =>  {len(D_members)} = {len(G_members)} - {len(G_next)}")
                else:
                    lines.append(f"G_{{{k+1}}} = ∅ (empty set)")
                lines.append("=" * 60)

                text1.delete(1.0, tk.END)
                text1.insert(tk.END, "\n".join(lines))
            except Exception as e:
                messagebox.showerror("Error", f"Error refreshing: {str(e)}")

        combo.bind('<<ComboboxSelected>>', lambda e: refresh())
        refresh()

        # ========== Tab 2: Duplicate Samples ==========
        tab2 = ttk.Frame(notebook)
        notebook.add(tab2, text="Duplicate Samples")

        dup_pairs = self.grading_result['duplicate_pairs_df']
        dup_groups = self.grading_result['duplicate_groups']

        if dup_pairs.empty:
            ttk.Label(tab2, text="✓ No duplicate samples detected (Mmmd=1)",
                      font=('Microsoft YaHei', 11, 'bold'), foreground='green').pack(pady=20)
        else:
            info_frame = ttk.Frame(tab2)
            info_frame.pack(fill=tk.X, padx=10, pady=5)
            ttk.Label(info_frame, text=f"Detected {len(dup_pairs)} duplicate sample pairs",
                      font=('Microsoft YaHei', 11, 'bold'), foreground='red').pack(side=tk.LEFT, padx=5)
            ttk.Label(info_frame, text=f"Duplicate groups: {len(dup_groups)} groups").pack(side=tk.LEFT, padx=5)

            text2 = scrolledtext.ScrolledText(tab2, wrap=tk.WORD, font=self.font_mono, width=80, height=36)
            text2.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

            lines = [
                "Duplicate Sample Detection Report (Mmmd = 1)",
                "=" * 60,
                f"Duplicate groups: {len(dup_groups)}",
                "=" * 60,
                "",
            ]

            for grp_idx, group in enumerate(dup_groups, 1):
                lines.append(f"[Duplicate Group {grp_idx}] Samples: {len(group)}")
                lines.append(f"  Members: {', '.join(group)}")
                lines.append("")

                group_pairs = []
                for _, row in dup_pairs.iterrows():
                    s1, s2 = row['SampleID_1'], row['SampleID_2']
                    if s1 in group and s2 in group:
                        group_pairs.append((s1, s2))

                if group_pairs:
                    lines.append(f"  Duplicate pairs ({len(group_pairs)} pairs):")
                    for s1, s2 in group_pairs:
                        lines.append(f"    {s1}  ↔  {s2}  (Mmmd=1)")
                lines.append("-" * 60)
                lines.append("")

            text2.insert(tk.END, "\n".join(lines))

        ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)

    def _view_duplicate_samples(self):
        """View duplicate samples"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        dup_pairs = self.grading_result['duplicate_pairs_df']
        dup_groups = self.grading_result['duplicate_groups']

        win = tk.Toplevel(self.root)
        win.title("Duplicate Sample Detection Report")
        win.geometry("650x720")
        win.transient(self.root)
        win.grab_set()

        if dup_pairs.empty:
            ttk.Label(win, text="✓ No duplicate samples detected (Mmmd=1)",
                      font=('Microsoft YaHei', 11, 'bold'), foreground='green').pack(pady=10)
            ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)
            return

        info_frame = ttk.Frame(win)
        info_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(info_frame, text=f"Detected {len(dup_pairs)} duplicate sample pairs",
                  font=('Microsoft YaHei', 11, 'bold'), foreground='red').pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Duplicate groups: {len(dup_groups)} groups").pack(side=tk.LEFT, padx=5)

        text = scrolledtext.ScrolledText(win, wrap=tk.WORD, font=self.font_mono, width=76, height=32)
        text.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        lines = [
            "Duplicate Sample Detection Report (Mmmd = 1)",
            "=" * 60,
            f"Duplicate groups: {len(dup_groups)}",
            "=" * 60,
            "",
        ]

        for grp_idx, group in enumerate(dup_groups, 1):
            lines.append(f"[Duplicate Group {grp_idx}] Samples: {len(group)}")
            lines.append(f"  Members: {', '.join(group)}")
            lines.append("")

            group_pairs = []
            for _, row in dup_pairs.iterrows():
                s1, s2 = row['SampleID_1'], row['SampleID_2']
                if s1 in group and s2 in group:
                    group_pairs.append((s1, s2))

            if group_pairs:
                lines.append(f"  Duplicate pairs ({len(group_pairs)} pairs):")
                for s1, s2 in group_pairs:
                    lines.append(f"    {s1}  ↔  {s2}  (Mmmd=1)")
            lines.append("-" * 60)
            lines.append("")

        text.insert(tk.END, "\n".join(lines))
        ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)

    def _view_d0_network(self):
        """View D0 network connectivity analysis"""
        if self.grading_result is None:
            messagebox.showwarning("Warning", "Please run grading analysis first!")
            return

        d0_members = self.grading_result['D_dict'].get(0, [])
        if not d0_members:
            messagebox.showinfo("Info", "D0 level is empty, no samples with Mmin=0.")
            return

        M_matrix = self.grading_result['M_matrix']
        sample_list = self.grading_result['sample_list']
        n_d0 = len(d0_members)

        d0_indices = [sample_list.index(m) for m in d0_members]

        # Union-find
        parent = list(range(n_d0))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        edges = []
        for i in range(n_d0):
            for j in range(i + 1, n_d0):
                val = M_matrix[d0_indices[i], d0_indices[j]]
                if not np.isnan(val) and val == 0:
                    union(i, j)
                    edges.append((d0_members[i], d0_members[j]))

        groups = {}
        for i in range(n_d0):
            root = find(i)
            groups.setdefault(root, []).append(i)

        components = sorted(groups.values(), key=lambda g: -len(g))

        # Calculate statistics
        comp_data = []
        for comp in components:
            comp_members = [d0_members[i] for i in comp]
            comp_data.append(comp_members)

        max_size = len(comp_data[0]) if comp_data else 0
        avg_size = sum(len(c) for c in comp_data) / len(comp_data) if comp_data else 0
        isolated = sum(1 for c in comp_data if len(c) == 1)

        win = tk.Toplevel(self.root)
        win.title("D0 Level Network Connectivity Analysis")
        win.geometry("700x720")
        win.transient(self.root)
        win.grab_set()

        # Statistics bar
        info_frame = ttk.LabelFrame(win, text="Statistics", padding=5)
        info_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(info_frame, text=f"D0 level total samples: {n_d0}").pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Connected edges: {len(edges)}").pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Connected components: {len(components)}").pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Max component: {max_size} items").pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Avg size: {avg_size:.1f} items").pack(side=tk.LEFT, padx=5)
        ttk.Label(info_frame, text=f"Isolated samples: {isolated}").pack(side=tk.LEFT, padx=5)

        # Left: connected component list
        paned = ttk.PanedWindow(win, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        left_frame = ttk.LabelFrame(paned, text="Connected Component List", padding=5)
        paned.add(left_frame, weight=1)

        cols = ("No.", "Sample Count", "Proportion")
        tree = ttk.Treeview(left_frame, columns=cols, show="headings", height=20)
        for c in cols:
            tree.heading(c, text=c)
        tree.column("No.", width=60, anchor='center')
        tree.column("Sample Count", width=70, anchor='center')
        tree.column("Proportion", width=70, anchor='center')
        tree.pack(fill=tk.BOTH, expand=True)

        # Right: member details
        right_frame = ttk.LabelFrame(paned, text="Component Member Details", padding=5)
        paned.add(right_frame, weight=2)

        detail_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD, font=self.font_mono, height=20)
        detail_text.pack(fill=tk.BOTH, expand=True)

        # Fill component list
        for comp_idx, comp_members in enumerate(comp_data, 1):
            pct = len(comp_members) / n_d0 * 100 if n_d0 > 0 else 0
            tree.insert("", tk.END, values=(comp_idx, len(comp_members), f"{pct:.1f}%"))

        def on_select(event):
            sel = tree.selection()
            if not sel:
                return
            idx = int(tree.item(sel[0], 'values')[0]) - 1
            comp_members = comp_data[idx]

            lines = [
                f"Connected Component {idx + 1}",
                "=" * 50,
                f"Sample count: {len(comp_members)}",
                "=" * 50,
                "",
            ]
            for m in comp_members:
                lines.append(f"  - {m}")
            lines.append("")
            lines.append("=" * 50)

            # Internal edges of this component
            comp_edges = []
            for e1, e2 in edges:
                if e1 in comp_members and e2 in comp_members:
                    comp_edges.append((e1, e2))
            if comp_edges:
                lines.append(f"Internal connected edges ({len(comp_edges)} edges):")
                for e1, e2 in comp_edges:
                    lines.append(f"  {e1}  ↔  {e2}")
            else:
                lines.append("Isolated sample, no connected edges")

            detail_text.delete(1.0, tk.END)
            detail_text.insert(tk.END, "\n".join(lines))

        tree.bind('<<TreeviewSelect>>', on_select)
        if comp_data:
            tree.selection_set(tree.get_children()[0])
            on_select(None)

        ttk.Button(win, text="Close", command=win.destroy).pack(pady=5)

    def _clear_grading_log(self):
        self.grading_log.configure(state=tk.NORMAL)
        self.grading_log.delete(1.0, tk.END)
        self.grading_log.configure(state=tk.DISABLED)

    # ==================== Common Methods ====================

    def _poll_queue(self):
        """Poll output queue in main thread and update log"""
        try:
            while True:
                text = self.output_queue.get_nowait()
                # Update currently visible tab log
                current_tab = self.notebook.index(self.notebook.select())
                if current_tab == 0:
                    log_widget = self.builder_log
                else:
                    log_widget = self.grading_log

                log_widget.configure(state=tk.NORMAL)
                log_widget.insert(tk.END, text)
                log_widget.see(tk.END)
                log_widget.configure(state=tk.DISABLED)
        except queue.Empty:
            pass
        self.root.after(100, self._poll_queue)

    def _on_close(self):
        self.root.destroy()
        sys.exit(0)


# =============================================================================
# Part 15: CLI and Program Entry
# =============================================================================

def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='Mendelian Mismatch Matrix Analysis System - v1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Launch GUI (auto-launched when no arguments)
  python MMMv1.py

  # Build matrix from raw genotype and auto-grade
  python MMMv1.py -i input.csv -o output_dir

  # Perform grading analysis on existing matrix only
  python MMMv1.py --grading -m MMM_Matrix.csv -o grading_output

  # Specify SSR marker
  python MMMv1.py -i input.csv -o output_dir -t SSR
        """
    )

    parser.add_argument('-i', '--input', help='Input raw genotype CSV file path')
    parser.add_argument('-o', '--output', help='Output directory path')
    parser.add_argument('-m', '--matrix', help='Existing MMM matrix file path (grading mode)')
    parser.add_argument('-t', '--type', choices=['SSR', 'SNP'], help='Marker type (optional, default auto-detect)')
    parser.add_argument('--no-grading', action='store_true', help='Do not auto-run grading analysis after matrix build')
    parser.add_argument('--grading', action='store_true', help='Execute grading analysis mode only')
    parser.add_argument('--no-threads', action='store_true', help='Disable multithreading')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')

    return parser


def main():
    """Main entry function"""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8")

    # Auto-launch GUI when no CLI arguments
    if len(sys.argv) <= 1:
        root = tk.Tk()
        app = MMMApp(root)
        root.mainloop()
        return 0

    parser = create_parser()
    args = parser.parse_args()

    # GUI mode
    if not args.input and not args.matrix:
        root = tk.Tk()
        app = MMMApp(root)
        root.mainloop()
        return 0

    # Grading analysis mode
    if args.grading or args.matrix:
        if not args.matrix or not args.output:
            parser.print_help()
            print("\nError: Grading mode requires --matrix and --output")
            return 2

        config = MMMConfig(
            input_file=args.matrix,
            output_dir=args.output,
            auto_grading=False
        )
        try:
            workflow = MMMWorkflow(config)
            result = workflow.run_grading_only(args.matrix)
            return 0
        except Exception as e:
            print(f"\n✗ Error: {str(e)}")
            import traceback
            traceback.print_exc()
            return 2

    # Matrix build mode
    if not args.input or not args.output:
        parser.print_help()
        print("\nError: Matrix build mode requires --input and --output")
        return 2

    config = MMMConfig(
        input_file=args.input,
        output_dir=args.output,
        marker_type=args.type,
        use_multithreading=not args.no_threads,
        auto_grading=not args.no_grading
    )

    try:
        workflow = MMMWorkflow(config)
        results = workflow.run()

        if results.validation_report.get("compliant", False):
            return 0
        else:
            print("\n⚠ Warning: Matrix did not pass all mathematical validations")
            return 1

    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return 2


if __name__ == "__main__":
    sys.exit(main())
