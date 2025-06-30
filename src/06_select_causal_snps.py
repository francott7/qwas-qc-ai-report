#!/usr/bin/env python3

"""
Select a random set of causal SNPs from a PLINK .raw genotype file for phenotype simulation.
Author: RFtt7

Inputs:
    --genotype-file: Path to the PLINK .raw file (e.g., data/processed/causal_genotypes.raw)
    --output-file:   Path to save the list of selected causal SNP IDs (e.g., data/processed/causal_snps.list)
    --n-causal-snps: Number of causal SNPs to select (default: 50)

Outputs:
    - A text file (output-file) containing one SNP ID per line, formatted for PLINK (e.g., rsID, without allele suffix)

Example usage:
    python src/06_select_causal_snps.py \
        --genotype-file data/processed/causal_genotypes.raw \
        --output-file data/processed/causal_snps.list \
        --n-causal-snps 50
"""

import numpy as np
import argparse
import logging
import subprocess
from utils import validate_file_exists

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def select_and_save_causal_snps(genotype_file: str, n_causal_snps: int, output_file: str):
    """
    Selects a random subset of SNPs to be 'causal' and saves their IDs to a file.
    """
    try:
        logger.info(f"Identifying SNP columns from {genotype_file}...")
        
        header_cmd = f"head -n 1 {genotype_file}"
        header_line = subprocess.check_output(header_cmd, shell=True).decode('utf-8').strip()
        all_cols = header_line.split()
        
        snp_cols = [col for col in all_cols if '_' in col and col not in ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        
        if len(snp_cols) < n_causal_snps:
            raise ValueError(f"Requested {n_causal_snps} causal SNPs, but only {len(snp_cols)} are available.")
            
        logger.info(f"Found {len(snp_cols)} SNPs. Selecting {n_causal_snps} to be causal.")
        
        np.random.seed(42)
        causal_snps = np.random.choice(snp_cols, n_causal_snps, replace=False)
        
        with open(output_file, 'w') as f:
            for snp_id in causal_snps:
                # The .raw file columns are named like 'SNP_ID_ALT', but PLINK needs just 'SNP_ID'
                base_snp_id = snp_id.rsplit('_', 1)[0]
                f.write(f"{base_snp_id}\n")
                
        logger.info(f"Saved {len(causal_snps)} causal SNP IDs to {output_file}")

    except Exception as e:
        logger.error(f"Failed to select and save causal SNPs: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Select a random set of causal SNPs for phenotype simulation.')
    parser.add_argument('--genotype-file', required=True, help='PLINK .raw file to read SNP IDs from.')
    parser.add_argument('--output-file', required=True, help='File to save the list of causal SNP IDs.')
    parser.add_argument('--n-causal-snps', type=int, default=50, help='Number of causal SNPs to select.')
    
    args = parser.parse_args()

    validate_file_exists(args.genotype_file, "Genotype file")
    
    select_and_save_causal_snps(
        genotype_file=args.genotype_file,
        n_causal_snps=args.n_causal_snps,
        output_file=args.output_file
    )
    
    return 0

if __name__ == "__main__":
    exit(main()) 