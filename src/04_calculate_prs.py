"""
Polygenic Risk Score Calculation Pipeline
Author: RFtt7
Copyright (c) 2025

Inputs:
    --vcf: Input VCF file with genotypes (e.g., data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz)
    --effect-sizes: GWAS effect sizes file (e.g., data/processed/gwas_effect_sizes_filtered_harmonized.txt)
    --metrics: (Optional) QC metrics file for PCA adjustment (e.g., results/01_qc/qc_metrics.json)
    --phenotype-file: (Optional) Phenotype file for regression (e.g., data/processed/simulated_phenotype.tsv)
    --output: Output file for PRS results (e.g., results/04_prs/prs_metrics.json)
    --standardize: (Optional) Standardize PRS to Z-scores

Outputs:
    - PRS results file (JSON or text, e.g., results/04_prs/prs_metrics.json)
    - (Optional) Regression results if phenotype file is provided

Example usage:
    python src/04_calculate_prs.py \
        --vcf data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --effect-sizes data/processed/gwas_effect_sizes_filtered_harmonized.txt \
        --metrics results/01_qc/qc_metrics.json \
        --output results/04_prs/prs_metrics.json \
        --standardize
"""

# Input: data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz (VCF), data/processed/gwas_effect_sizes_filtered_harmonized.txt (GWAS summary stats), results/01_qc/qc_metrics.json (QC metrics)

#!/usr/bin/env python3

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
import logging
import subprocess
from typing import Dict, List, Tuple, Optional, Any
from utils import validate_file_exists, safe_json_dump

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def gt_to_dosage(gt: str) -> float:
    """Convert genotype string to dosage value."""
    if gt in {'0/0', '0|0'}:
        return 0.0
    elif gt in {'0/1', '1/0', '0|1', '1|0'}:
        return 1.0
    elif gt in {'1/1', '1|1'}:
        return 2.0
    else:
        return np.nan

def create_position_file(effects: pd.DataFrame, output_file: str = 'snps.pos') -> None:
    """Create a position file for bcftools in chr\tstart\tend format."""
    # Create a DataFrame with chr, start, end columns
    pos_df = pd.DataFrame({
        'chr': effects['CHR'].astype(str),
        'start': effects['BP'].astype(str),
        'end': effects['BP'].astype(str)
    })
    # Save to file in tab-separated format without header
    pos_df.to_csv(output_file, sep='\t', header=False, index=False)
    logger.info(f"Created position file with {len(pos_df)} positions")

def stream_prs_from_vcf(vcf_file: str, effects: pd.DataFrame) -> Tuple[np.ndarray, List[str]]:
    """Stream PRS calculation from VCF file to prevent memory issues.
    
    Args:
        vcf_file: Path to VCF file
        effects: DataFrame with effect sizes
        
    Returns:
        Tuple of (PRS array, sample IDs)
    """
    beta_map = dict(zip(effects['pos_key'], effects['BETA']))
    sample_cmd = f"bcftools query -l {vcf_file}"
    samples = subprocess.check_output(sample_cmd, shell=True, text=True).strip().split('\n')
    prs = np.zeros(len(samples))
    used = 0

    # Create position file for bcftools
    create_position_file(effects)
    
    # Stream variant by variant
    fmt = '%CHROM\t%POS[\t%GT]\n'
    query_cmd = ["bcftools", "query", "-f", fmt, "-R", "snps.pos", vcf_file]
    
    try:
        with subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                fields = line.rstrip().split('\t')
                pos_key = f"{fields[0]}:{fields[1]}"
                beta = beta_map.get(pos_key)
                if beta is None:  # safety check
                    continue
                dos = np.array([gt_to_dosage(gt) for gt in fields[2:]])
                valid_dos = dos[~np.isnan(dos)]
                if len(valid_dos) > 0:  # only impute if we have valid dosages
                    dos[np.isnan(dos)] = valid_dos.mean()
                    prs += dos * beta
                    used += 1
                
                if used % 1000 == 0:  # Progress logging
                    logger.info(f"Processed {used} SNPs")
                    
        logger.info(f"Streaming PRS calculation completed using {used} SNPs")
        return prs, samples
        
    except Exception as e:
        logger.error(f"Error in streaming PRS calculation: {e}")
        return np.array([]), []
    finally:
        # Cleanup
        if Path('snps.pos').exists():
            Path('snps.pos').unlink()

def load_gwas_effect_sizes(effect_sizes_file: str) -> Optional[pd.DataFrame]:
    """Load GWAS effect sizes from file.
    
    Expected format: SNP, CHR, BP, A1, A2, BETA, SE, P
    """
    if not validate_file_exists(effect_sizes_file, "GWAS effect sizes"):
        return None
        
    try:
        df = pd.read_csv(effect_sizes_file, sep='\t')
        required_columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            logger.error(f"Missing required columns in effect sizes file: {missing_columns}")
            return None
            
        # Add position key for matching
        df['pos_key'] = df['CHR'].astype(str) + ':' + df['BP'].astype(str)
        logger.info(f"Loaded {len(df)} effect sizes")
        return df
    except Exception as e:
        logger.error(f"Error loading GWAS effect sizes: {e}")
        return None

def standardize_prs(prs_values: np.ndarray) -> Optional[np.ndarray]:
    """Standardize PRS values to Z-scores."""
    try:
        if len(prs_values) == 0:
            logger.error("No PRS values to standardize")
            return None
            
        mean_prs = np.mean(prs_values)
        std_prs = np.std(prs_values)
        
        if std_prs == 0:
            logger.error("Cannot standardize PRS: standard deviation is zero")
            return None
            
        standardized = (prs_values - mean_prs) / std_prs
        logger.info("PRS values standardized successfully")
        return standardized
    except Exception as e:
        logger.error(f"Error standardizing PRS values: {e}")
        return None

def load_pca_from_metrics(metrics_file: str, n_pcs: int = 10) -> Optional[Dict[str, Any]]:
    """Load PCA results (explained variance and PCs) from QC metrics JSON."""
    if not validate_file_exists(metrics_file, "QC metrics"):
        return None
    try:
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)
        pca = metrics.get('pca', {})
        if not pca or 'pcs' not in pca:
            logger.warning("No PCA results found in metrics file.")
            return None
        return {
            'explained_variance': pca.get('explained_variance', []),
            'pcs': pca.get('pcs', [])[:n_pcs]
        }
    except Exception as e:
        logger.error(f"Error loading PCA from metrics: {e}")
        return None

def run_prs_regression(phenotype_file: str, prs_values: np.ndarray, pcs: np.ndarray) -> Optional[Dict[str, Any]]:
    """Run regression of phenotype ~ PRS + PCs if phenotype data is available."""
    try:
        if not validate_file_exists(phenotype_file, "Phenotype"):
            logger.warning("Phenotype file not found. Skipping regression.")
            return None
        pheno_df = pd.read_csv(phenotype_file, sep=None, engine='python')
        if 'phenotype' not in pheno_df.columns:
            logger.error("Phenotype column not found in phenotype file.")
            return None
        from sklearn.linear_model import LinearRegression
        X = np.column_stack([prs_values, pcs])
        y = pheno_df['phenotype'].values
        model = LinearRegression().fit(X, y)
        r2 = model.score(X, y)
        return {'r2': r2, 'coef': model.coef_.tolist(), 'intercept': model.intercept_}
    except Exception as e:
        logger.error(f"Error in PRS regression: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Calculate PRS from genotype data and GWAS effect sizes')
    parser.add_argument('--vcf', required=True, help='Input VCF file with genotypes')
    parser.add_argument('--effect-sizes', required=True, help='Input file with GWAS effect sizes')
    parser.add_argument('--output', required=True, help='Output file for PRS results')
    parser.add_argument('--metrics', help='Input file with QC metrics')
    parser.add_argument('--standardize', action='store_true',
                       help='Standardize PRS to Z-scores')
    parser.add_argument('--phenotype-file', help='Input file with phenotype data')
    
    args = parser.parse_args()
    
    # Load effect sizes first (smaller file)
    effect_sizes_df = load_gwas_effect_sizes(args.effect_sizes)
    if effect_sizes_df is None:
        logger.error("Failed to load effect sizes. Exiting.")
        return
    
    # Calculate PRS using streaming approach
    prs, sample_ids = stream_prs_from_vcf(args.vcf, effect_sizes_df)
    if len(prs) == 0:
        logger.error("PRS calculation failed. Exiting.")
        return
        
    # Standardize if requested
    standardized_prs = standardize_prs(prs) if args.standardize else None
    
    # Load PCA from metrics if available
    pca_results = load_pca_from_metrics(args.metrics, n_pcs=10) if args.metrics else None
    pcs = np.array(pca_results['pcs']) if pca_results else np.empty((0, 10))
    
    # Run regression if phenotype data is available
    regression_results = run_prs_regression(args.phenotype_file, prs, pcs) if args.phenotype_file else {}
    
    # Prepare and save results
    results = {
        'raw_prs': prs.tolist(),
        'standardized_prs': standardized_prs.tolist() if standardized_prs is not None else None,
        'sample_ids': sample_ids,
        'n_snps_used': len(effect_sizes_df),
        'pcs': pcs.tolist() if pcs.size else [],
        'regression': regression_results or {}
    }
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save results
    safe_json_dump(results, str(output_path))
    logger.info(f"Results saved to {args.output}")

if __name__ == '__main__':
    main() 