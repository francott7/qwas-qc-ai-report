#!/usr/bin/env python3

"""
Simulate a case-control phenotype for PRS validation using real genotype data.
Author: RFtt7

Inputs:
    --genotype-file: PLINK .raw file with genotype dosages (e.g., data/processed/causal_genotypes.raw)
    --output-file: Path to save the simulated phenotype file (e.g., data/processed/simulated_phenotype.tsv)
    --heritability: Heritability of the simulated trait (default: 0.5)
    --prevalence: Case prevalence in the population (default: 0.2)

Outputs:
    - Simulated phenotype file (TSV) with columns: FID, IID, PHENO

Example usage:
    python src/07_simulate_phenotype.py \
        --genotype-file data/processed/causal_genotypes.raw \
        --output-file data/processed/simulated_phenotype.tsv \
        --heritability 0.5 \
        --prevalence 0.2
"""

import pandas as pd
import numpy as np
import argparse
import os
import logging
from utils import validate_file_exists, safe_save_df

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def simulate_phenotype(genotype_file: str, n_causal_snps: int, heritability: float, prevalence: float) -> pd.DataFrame:
    """
    Simulates a case-control phenotype based on real genotype data.
    """
    try:
        logger.info(f"Loading genotype data from {genotype_file}")
        geno_df = pd.read_csv(genotype_file, delim_whitespace=True)

        snp_cols = [col for col in geno_df.columns if '_' in col and col not in ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        
        if not snp_cols:
            raise ValueError("No SNP columns found in the genotype file.")
        
        # All SNPs in this file are causal
        causal_snps = snp_cols
        n_causal_snps = len(causal_snps)
        logger.info(f"Using {n_causal_snps} causal SNPs from the file.")

        # Assign effect sizes (betas)
        np.random.seed(42)
        effect_sizes = np.random.normal(0, 1, n_causal_snps)

        # Calculate the genetic component (true PRS)
        dosages = geno_df[causal_snps].values
        # Handle potential missing values that PLINK might output as NA
        dosages = np.nan_to_num(dosages)
        genetic_component = np.dot(dosages, effect_sizes)

        # Calculate environmental component
        genetic_variance = np.var(genetic_component)
        environmental_variance = genetic_variance * (1 / heritability - 1) if heritability > 0 else 1
        environmental_noise = np.random.normal(0, np.sqrt(environmental_variance), len(geno_df))

        # Create the full liability score
        liability_score = genetic_component + environmental_noise

        # Threshold to create a binary case/control phenotype
        threshold = np.quantile(liability_score, 1 - prevalence)
        phenotype = (liability_score >= threshold).astype(int)

        # Create the final phenotype DataFrame
        pheno_df = geno_df[['FID', 'IID']].copy()
        pheno_df['PHENO'] = phenotype
        
        logger.info(f"Simulation complete. Generated {pheno_df['PHENO'].sum()} cases and {len(pheno_df) - pheno_df['PHENO'].sum()} controls.")
        
        actual_h2 = genetic_variance / np.var(liability_score) if np.var(liability_score) > 0 else 0
        logger.info(f"Target heritability: {heritability:.3f}, Actual simulated heritability: {actual_h2:.3f}")

        return pheno_df

    except Exception as e:
        logger.error(f"Phenotype simulation failed: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Simulate a case-control phenotype for PRS validation.')
    parser.add_argument('--genotype-file', required=True, help='PLINK .raw file with genotype dosages.')
    parser.add_argument('--output-file', required=True, help='Path to save the simulated phenotype file.')
    parser.add_argument('--heritability', type=float, default=0.5, help='Heritability of the simulated trait (0 to 1).')
    parser.add_argument('--prevalence', type=float, default=0.2, help='Case prevalence in the population (0 to 1).')
    
    args = parser.parse_args()
    
    validate_file_exists(args.genotype_file, "Genotype file")

    # In this new workflow, n_causal_snps is determined by the input file
    # We pass a placeholder value, it will be ignored.
    simulate_phenotype(
        genotype_file=args.genotype_file,
        n_causal_snps=-1, # Placeholder
        heritability=args.heritability,
        prevalence=args.prevalence
    )

    simulated_pheno_df = simulate_phenotype(
        genotype_file=args.genotype_file,
        n_causal_snps=-1,  # Placeholder, will be determined from file
        heritability=args.heritability,
        prevalence=args.prevalence
    )
    
    safe_save_df(simulated_pheno_df, args.output_file, index=False, sep='\t')
    logger.info(f"Simulated phenotype data saved to {args.output_file}")
    return 0

if __name__ == "__main__":
    exit(main()) 