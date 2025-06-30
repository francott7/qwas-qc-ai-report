"""
GWAS Quality Control Pipeline using Hail
Author: RFtt7
Copyright (c) 2025

Inputs:
    --vcf: Path to input VCF file (e.g., data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz)

Outputs:
    - QC metrics JSON file (e.g., results/01_qc/qc_metrics.json)
    - PCA results (in JSON)
    - Summary plots (optional, e.g., in results/01_qc/)

Example usage:
    python src/01_qc_hail.py \
        --vcf data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --output results/01_qc/qc_metrics.json
"""

# Input: data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

#!/usr/bin/env python3

import argparse
import json
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, Any, Optional, Tuple
from utils import (
    check_dependency,
    validate_file_exists,
    validate_vcf_format,
    safe_json_dump,
    get_file_size,
    check_disk_space
)
import tempfile
import os
import multiprocessing
import random
from sklearn.decomposition import PCA
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_bcftools_stats(vcf_file: str) -> Optional[str]:
    """Run bcftools stats on the VCF file."""
    if not check_dependency('bcftools'):
        return None
        
    cmd = f"bcftools stats {vcf_file}"
    logger.info(f"Running command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        result.check_returncode()
        logger.info("bcftools stats completed successfully")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running bcftools stats: {e.stderr}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error running bcftools stats: {e}")
        return None

def parse_bcftools_stats(stats_output: str) -> Dict[str, Any]:
    """Parse bcftools stats output into a dictionary."""
    if not stats_output:
        logger.error("No stats output received from bcftools")
        return {}
        
    stats = {}
    current_section = None
    
    try:
        for line in stats_output.split('\n'):
            if line.startswith('#'):
                continue
            if line.startswith('SN'):
                fields = line.split('\t')
                if len(fields) >= 3:
                    key = fields[2].strip().rstrip(':')
                    value = fields[3].strip()
                    stats[key] = value
                    logger.debug(f"Found metric: {key} = {value}")
            elif line.startswith('TSTV'):
                logger.debug(f"Processing TSTV line: {line}")
                fields = [f for f in line.split() if f]
                
                if len(fields) >= 8:
                    try:
                        stats['ts_count'] = int(fields[2])
                        stats['tv_count'] = int(fields[3])
                        stats['ts_tv_ratio'] = float(fields[4])
                        stats['ts_count_1st_alt'] = int(fields[5])
                        stats['tv_count_1st_alt'] = int(fields[6])
                        stats['ts_tv_ratio_1st_alt'] = float(fields[7])
                        logger.info(f"TSTV stats parsed successfully")
                    except (ValueError, IndexError) as e:
                        logger.error(f"Error parsing TSTV line: {e}")
                        logger.error(f"Problematic line: {line}")
                        logger.error(f"Fields: {fields}")
    except Exception as e:
        logger.error(f"Error parsing bcftools stats: {e}")
        return {}
    
    return stats

def calculate_dosage_metrics(vcf_file: str) -> Dict[str, Any]:
    """Calculate dosage-specific metrics from VCF file."""
    try:
        # Run bcftools query to get dosage information
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/DS\\n' {vcf_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        result.check_returncode()
        
        # Parse dosage data
        dosages = []
        for line in result.stdout.split('\n'):
            if line.strip():
                fields = line.split('\t')
                if len(fields) >= 6:
                    try:
                        dosage = float(fields[5])
                        dosages.append(dosage)
                    except ValueError:
                        continue
        
        if not dosages:
            return {
                'avg_quality': 0.0,
                'snps_below_threshold': 0,
                'distribution': {'0': 0, '1': 0, '2': 0}
            }
        
        # Calculate metrics
        dosages = np.array(dosages)
        return {
            'avg_quality': float(np.mean(dosages)),
            'snps_below_threshold': int(np.sum(dosages < 0.1)),
            'distribution': {
                '0': int(np.sum(dosages < 0.5)),
                '1': int(np.sum((dosages >= 0.5) & (dosages < 1.5))),
                '2': int(np.sum(dosages >= 1.5))
            }
        }
    except Exception as e:
        logger.warning(f"Dosage metrics unavailable: {e}")
        return {
            'avg_quality': 0.0,
            'snps_below_threshold': 0,
            'distribution': {'0': 0, '1': 0, '2': 0}
        }

def process_gt_chunk(lines):
    total = 0
    non_missing = 0
    het = 0
    non_ref = 0
    alt_alleles = 0
    ref_alleles = 0
    for gt in lines:
        gt = gt.strip()
        total += 1
        if gt == './.':
            continue
        non_missing += 1
        # Heterozygous: 0/1, 1/0, 0|1, 1|0
        if gt in {'0/1', '1/0', '0|1', '1|0'}:
            het += 1
        # Non-ref: not 0/0 or 0|0
        if gt not in {'0/0', '0|0'}:
            non_ref += 1
        # Count alleles
        if '/' in gt or '|' in gt:
            sep = '/' if '/' in gt else '|'
            alleles = gt.split(sep)
            for a in alleles:
                if a == '.':
                    continue
                elif a == '0':
                    ref_alleles += 1
                elif a == '1':
                    alt_alleles += 1
                # ignore multi-allelic for now
    return total, non_missing, het, non_ref, alt_alleles, ref_alleles

def parallel_gt_metrics(gt_tmp_path, n_workers=8, chunk_size=1000000):
    pool = multiprocessing.Pool(n_workers)
    results = []
    with open(gt_tmp_path) as f:
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                results.append(pool.apply_async(process_gt_chunk, (chunk,)))
                chunk = []
        if chunk:
            results.append(pool.apply_async(process_gt_chunk, (chunk,)))
    pool.close()
    pool.join()
    # Aggregate
    total = non_missing = het = non_ref = alt_alleles = ref_alleles = 0
    for r in results:
        t, nm, h, nr, alt, ref = r.get()
        total += t
        non_missing += nm
        het += h
        non_ref += nr
        alt_alleles += alt
        ref_alleles += ref
    return total, non_missing, het, non_ref, alt_alleles, ref_alleles

def run_pca_on_vcf(vcf_file, n_snps=10000, n_components=10, random_seed=42):
    """Extract a random subset of SNPs and run PCA using scikit-learn."""
    import tempfile
    import subprocess
    import pandas as pd
    random.seed(random_seed)
    try:
        # 1. Get all variant positions
        cmd_sites = f"bcftools query -f '%CHROM\t%POS\n' {vcf_file}"
        result = subprocess.run(cmd_sites, shell=True, capture_output=True, text=True)
        result.check_returncode()
        sites = [tuple(line.split('\t')) for line in result.stdout.strip().split('\n') if line]
        if len(sites) <= n_snps:
            selected_sites = sites
        else:
            selected_sites = random.sample(sites, n_snps)
        # 2. Extract GTs for selected sites
        site_str = ' '.join([f"{chrom}:{pos}" for chrom, pos in selected_sites])
        with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.gt') as gt_tmp:
            gt_tmp_path = gt_tmp.name
            cmd_query = f"bcftools query -r {','.join([f'{chrom}:{pos}' for chrom, pos in selected_sites])} -f '[%GT\t]' {vcf_file} > {gt_tmp_path}"
            subprocess.run(cmd_query, shell=True, check=True)
        # 3. Read GTs and convert to numeric
        with open(gt_tmp_path) as f:
            lines = [line.strip().split('\t') for line in f if line.strip()]
        os.remove(gt_tmp_path)
        if not lines:
            return None, None, "No GTs found for PCA."
        # Transpose: rows=samples, cols=variants
        gt_matrix = np.array(lines, dtype=object).T
        # Convert GTs to numeric (0,1,2), missing as np.nan
        def gt_to_num(gt):
            if gt in {'0/0', '0|0'}:
                return 0
            elif gt in {'0/1', '1/0', '0|1', '1|0'}:
                return 1
            elif gt in {'1/1', '1|1'}:
                return 2
            else:
                return np.nan
        gt_numeric = np.vectorize(gt_to_num)(gt_matrix)
        # Impute missing with mean for each SNP
        from sklearn.impute import SimpleImputer
        imputer = SimpleImputer(strategy='mean')
        gt_numeric = imputer.fit_transform(gt_numeric)
        # 4. Run PCA
        pca = PCA(n_components=min(n_components, gt_numeric.shape[1], gt_numeric.shape[0]))
        pcs = pca.fit_transform(gt_numeric)
        explained = pca.explained_variance_ratio_.tolist()
        pcs_list = pcs[:, :n_components].tolist()
        return explained, pcs_list, None
    except Exception as e:
        return None, None, str(e)

def calculate_qc_metrics(vcf_file: str, output_file: str) -> bool:
    """Calculate QC metrics using bcftools."""
    # Validate input file
    if not validate_file_exists(vcf_file, "VCF"):
        logger.critical("Input VCF file does not exist or is not accessible.")
        safe_json_dump({"error": "Input VCF file does not exist or is not accessible."}, output_file)
        return False
    
    if not validate_vcf_format(vcf_file):
        logger.critical("Input VCF file format is invalid.")
        safe_json_dump({"error": "Input VCF file format is invalid."}, output_file)
        return False
    
    # Check disk space (estimate 2x input file size for processing)
    file_size = get_file_size(vcf_file)
    if file_size and not check_disk_space(str(Path(output_file).parent), file_size * 2 / (1024 * 1024)):
        logger.critical("Insufficient disk space for processing.")
        safe_json_dump({"error": "Insufficient disk space for processing."}, output_file)
        return False
    
    # Run bcftools stats
    stats_output = run_bcftools_stats(vcf_file)
    if not stats_output:
        logger.critical("Failed to run bcftools stats.")
        safe_json_dump({"error": "Failed to run bcftools stats."}, output_file)
        return False
    
    # Parse stats
    stats = parse_bcftools_stats(stats_output)
    if not stats:
        logger.critical("Failed to parse bcftools stats output.")
        safe_json_dump({"error": "Failed to parse bcftools stats output."}, output_file)
        return False
    
    # Log and check key counts
    n_variants = int(stats.get('number of records', '0'))
    n_samples = int(stats.get('number of samples', '0'))
    n_snps = int(stats.get('number of SNPs', '0'))
    logger.info(f"Variants: {n_variants}, Samples: {n_samples}, SNPs: {n_snps}")
    qc_summary = {
        'variants': n_variants,
        'samples': n_samples,
        'snps': n_snps,
        'warnings': [],
        'errors': []
    }
    if n_variants == 0 or n_samples == 0 or n_snps == 0:
        msg = "Zero variants, samples, or SNPs detected. Check VCF integrity and pipeline import steps."
        logger.critical(msg)
        qc_summary['errors'].append(msg)
        safe_json_dump({'qc_summary': qc_summary}, output_file)
        return False
    
    # --- Optimized: Run bcftools query once, process temp file for all GT-based metrics ---
    try:
        with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.gt') as gt_tmp:
            gt_tmp_path = gt_tmp.name
            cmd_query = f"bcftools query -f '[%GT\\n]' {vcf_file} > {gt_tmp_path}"
            subprocess.run(cmd_query, shell=True, check=True)
        # Parallelized metrics
        total_gts, non_missing_gts, n_het, n_non_ref, alt_alleles, ref_alleles = parallel_gt_metrics(gt_tmp_path)
        call_rate = non_missing_gts / total_gts if total_gts else 0
        het_rate = n_het / non_missing_gts if non_missing_gts else 0
        non_ref_rate = n_non_ref / non_missing_gts if non_missing_gts else 0
        allele_freq = alt_alleles / (alt_alleles + ref_alleles) if (alt_alleles + ref_alleles) else 0
        logger.info(f"[Parallel] GTs: {total_gts}, Non-missing: {non_missing_gts}, Het: {n_het}, Non-ref: {n_non_ref}, Alt alleles: {alt_alleles}, Ref alleles: {ref_alleles}")
        # Clean up temp file
        os.remove(gt_tmp_path)
    except Exception as e:
        logger.error(f"Error in parallel GT metrics: {e}")
        call_rate = het_rate = non_ref_rate = allele_freq = 0
        n_het = n_non_ref = alt_alleles = ref_alleles = 0
        qc_summary['warnings'].append(f"Parallel GT metrics calculation failed: {e}")
    finally:
        try:
            if 'gt_tmp_path' in locals() and os.path.exists(gt_tmp_path):
                os.remove(gt_tmp_path)
        except Exception as cleanup_e:
            logger.warning(f"Could not remove temp GT file: {cleanup_e}")
    
    # Calculate dosage metrics
    dosage_metrics = calculate_dosage_metrics(vcf_file)
    
    # Compile all metrics
    metrics = {
        'number_of_variants': n_variants,
        'number_of_samples': n_samples,
        'ts_tv_ratio': float(stats.get('ts_tv_ratio', '0')),
        'ts_count': int(stats.get('ts_count', '0')),
        'tv_count': int(stats.get('tv_count', '0')),
        'ts_tv_ratio_1st_alt': float(stats.get('ts_tv_ratio_1st_alt', '0')),
        'ts_count_1st_alt': int(stats.get('ts_count_1st_alt', '0')),
        'tv_count_1st_alt': int(stats.get('tv_count_1st_alt', '0')),
        'number_of_snps': n_snps,
        'number_of_indels': int(stats.get('number of indels', '0')),
        'number_of_multiallelic_sites': int(stats.get('number of multiallelic sites', '0')),
        'number_of_multiallelic_SNP_sites': int(stats.get('number of multiallelic SNP sites', '0')),
        'number_of_multiallelic_indel_sites': int(stats.get('number of multiallelic indel sites', '0')),
        'sample_qc': {
            'call_rate': call_rate,
            'n_het': n_het,
            'n_non_ref': n_non_ref
        },
        'variant_qc': {
            'call_rate': call_rate,
            'allele_frequency': allele_freq
        },
        'dosage_qc': dosage_metrics,
        'qc_summary': qc_summary
    }
    
    # Save metrics to JSON
    if not safe_json_dump(metrics, output_file):
        logger.critical("Failed to write QC metrics to output JSON.")
        return False
    logger.info("QC metrics calculation completed successfully.")

    logger.info("Starting PCA calculation on genotype matrix...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        explained, pcs, pca_err = run_pca_on_vcf(vcf_file)
    pca_result = {}
    if explained is not None:
        logger.info(f"PCA explained variance: {explained}")
        pca_result = {'explained_variance': explained, 'pcs': pcs}
    else:
        logger.warning(f"PCA calculation failed: {pca_err}")
        pca_result = {'explained_variance': [], 'pcs': [], 'warning': pca_err}

    # Prepare output dict (ensure all metrics are included)
    output = {
        'number_of_variants': n_variants,
        'number_of_samples': n_samples,
        'ts_tv_ratio': float(stats.get('ts_tv_ratio', '0')),
        'ts_count': int(stats.get('ts_count', '0')),
        'tv_count': int(stats.get('tv_count', '0')),
        'ts_tv_ratio_1st_alt': float(stats.get('ts_tv_ratio_1st_alt', '0')),
        'ts_count_1st_alt': int(stats.get('ts_count_1st_alt', '0')),
        'tv_count_1st_alt': int(stats.get('tv_count_1st_alt', '0')),
        'number_of_snps': n_snps,
        'number_of_indels': int(stats.get('number of indels', '0')),
        'number_of_multiallelic_sites': int(stats.get('number of multiallelic sites', '0')),
        'number_of_multiallelic_SNP_sites': int(stats.get('number of multiallelic SNP sites', '0')),
        'number_of_multiallelic_indel_sites': int(stats.get('number of multiallelic indel sites', '0')),
        'sample_qc': {
            'call_rate': call_rate,
            'n_het': n_het,
            'n_non_ref': n_non_ref
        },
        'variant_qc': {
            'call_rate': call_rate,
            'allele_frequency': allele_freq
        },
        'dosage_qc': dosage_metrics,
        'qc_summary': qc_summary,
        'pca': pca_result
    }
    safe_json_dump(output, output_file)
    logger.info("QC metrics (including PCA) written to output JSON.")
    return True

def main():
    parser = argparse.ArgumentParser(description='Calculate QC metrics for a VCF file using bcftools')
    parser.add_argument('--input', required=True, help='Input VCF file')
    parser.add_argument('--file-type', choices=['vcf', 'bgen'], required=True, help='Input file type')
    parser.add_argument('--output', required=True, help='Output JSON file for QC metrics')
    
    args = parser.parse_args()
    
    if args.file_type == 'vcf':
        success = calculate_qc_metrics(args.input, args.output)
        if not success:
            logger.error("Failed to calculate QC metrics")
            exit(1)
    else:
        logger.error("BGEN file support not implemented yet")
        exit(1)

if __name__ == '__main__':
    main() 