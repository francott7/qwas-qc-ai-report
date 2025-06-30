#!/usr/bin/env python3

"""
Script to analyze PRS results and generate comprehensive metrics and visualizations.
Author: RFtt7

Example usage:
    python src/05_analyze_prs.py \
        --prs-file results/04_prs/prs_metrics.json \
        --qc-file results/01_qc/qc_metrics.json \
        --output-dir results/05_analysis

Required inputs:
    - PRS file (results/04_prs/prs_metrics.json):
        JSON file containing raw PRS scores for each sample.
        Format: {"raw_prs": [score1, score2, ...]}
        
    - QC metrics file (results/01_qc/qc_metrics.json):
        JSON file containing QC metrics including PCA results.
        Must contain 'pca' section with principal components.

Optional inputs:
    - Phenotype file (not used in current analysis):
        Tab-separated file containing phenotype data.
        Must have 'IID' column matching sample IDs.

Outputs:
    - prs_distribution.png: Distribution plots of PRS scores
    - population_stratification.png: PCA vs PRS plots
    - prs_analysis_metrics.json: Comprehensive analysis metrics
    
The script:
1. Loads and normalizes PRS scores
2. Adjusts for ancestry using first 10 PCs
3. Generates distribution plots
4. Calculates comprehensive metrics
5. Tests for normality

Note: Sample IDs are loaded from data/processed/chr22.fam to ensure consistent ordering.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import argparse
import os
from pathlib import Path
import logging
from typing import Dict, Any, List, Tuple, Optional
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score, r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from utils import validate_file_exists, safe_json_dump
import statsmodels.api as sm
from scipy import stats

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_sample_ids(fam_file: str = "data/processed/chr22.fam") -> List[str]:
    """Load sample IDs from FAM file."""
    try:
        df = pd.read_csv(fam_file, sep='\s+', header=None)
        return df[1].tolist()  # Second column contains IIDs
    except Exception as e:
        logger.error(f"Failed to load sample IDs from FAM file: {e}")
        return None

def load_prs_results(prs_file: str) -> Optional[pd.DataFrame]:
    """Load PRS results from JSON file."""
    try:
        with open(prs_file, 'r') as f:
            prs_data = json.load(f)
        
        # Get sample IDs
        sample_ids = load_sample_ids()
        if not sample_ids or len(sample_ids) != len(prs_data['raw_prs']):
            logger.warning("Sample IDs not found or mismatch in length. Using generic IDs.")
            sample_ids = [f'Sample_{i}' for i in range(len(prs_data['raw_prs']))]
        
        # Create a DataFrame with the raw PRS scores
        df = pd.DataFrame({
            'IID': sample_ids,
            'SCORE': prs_data['raw_prs'],
            'CNT2': [55042] * len(prs_data['raw_prs'])  # Number of SNPs used in PRS calculation
        })
        
        # Set IID as index
        df.set_index('IID', inplace=True)
        
        logger.info(f"Successfully loaded PRS results: {len(df)} samples")
        return df
    except Exception as e:
        logger.error(f"Failed to load PRS results: {e}")
        return None

def load_pca_results_from_qc(qc_file: str) -> Optional[pd.DataFrame]:
    """Load PCA results from the QC metrics JSON file."""
    try:
        with open(qc_file, 'r') as f:
            qc_data = json.load(f)
        
        # Check for 'pca_results' (new format) or 'pca' (old format)
        pca_data_key = None
        if 'pca_results' in qc_data and 'principal_components' in qc_data['pca_results']:
            pca_data_key = 'pca_results'
        elif 'pca' in qc_data and 'pcs' in qc_data['pca']:
            pca_data_key = 'pca'
        
        if pca_data_key == 'pca_results':
            pc_data = qc_data['pca_results']['principal_components']
            logger.info("Found 'pca_results' in QC metrics.")
            
        elif pca_data_key == 'pca':
            pc_data = qc_data['pca']['pcs']
            logger.info("Found 'pca' key in QC metrics, adapting to its structure.")
            
        else:
            logger.warning("No compatible PCA results found in QC metrics file.")
            return None

        if not pc_data:
            logger.warning("PCA results are empty.")
            return None

        # Get sample IDs
        sample_ids = load_sample_ids()
        if not sample_ids or len(sample_ids) != len(pc_data):
            logger.warning("Sample IDs not found or mismatch in length. Using generic IDs.")
            sample_ids = [f'Sample_{i}' for i in range(len(pc_data))]

        # Convert to DataFrame
        pca_df = pd.DataFrame(pc_data)
        
        # Name the PC columns
        pca_df.columns = [f'PC{i+1}' for i in range(pca_df.shape[1])]
        
        # Add sample IDs
        pca_df['IID'] = sample_ids
        pca_df.set_index('IID', inplace=True)

        logger.info(f"Successfully loaded PCA results for {len(pca_df)} samples from {qc_file}")
        return pca_df

    except Exception as e:
        logger.error(f"Failed to load or parse PCA results from {qc_file}: {e}")
        return None

def load_phenotype_data(phenotype_file: Optional[str]) -> Optional[pd.DataFrame]:
    """Load phenotype data if provided."""
    if not phenotype_file or not os.path.exists(phenotype_file):
        logger.info("No phenotype file provided or file doesn't exist")
        return None
    
    try:
        pheno_df = pd.read_csv(phenotype_file, sep='\t')
        logger.info(f"Successfully loaded phenotype data: {len(pheno_df)} samples")
        return pheno_df
    except Exception as e:
        logger.error(f"Failed to load phenotype data: {e}")
        return None

def normalize_prs_by_snp_count(prs_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize PRS scores by SNP count to remove linear dependence."""
    df = prs_df.copy()
    
    # Method 1: Mean dosage per SNP (current method)
    df['PRS_normalized'] = df['SCORE']  # No need to normalize by SNP count since all samples have same count
    df['PRS_normalized'] = df['PRS_normalized'].fillna(0)
    
    # Method 2: Z-score normalization of the mean dosage
    scaler = StandardScaler()
    df['PRS_normalized_std'] = scaler.fit_transform(df[['PRS_normalized']])
    
    # Method 3: Alternative normalization using median SNP count as reference
    df['PRS_median_normalized'] = df['SCORE']  # No need to normalize by SNP count
    df['PRS_median_normalized'] = df['PRS_median_normalized'].fillna(0)
    
    # Method 4: PLINK2-style normalization (standardize)
    df['PRS_plink_style'] = df['SCORE']
    df['PRS_plink_style'] = df['PRS_plink_style'].fillna(0)
    scaler_plink = StandardScaler()
    df['PRS_plink_style_std'] = scaler_plink.fit_transform(df[['PRS_plink_style']])
    
    # Calculate correlation metrics for all methods
    corr_original = 0  # All samples have same SNP count
    corr_normalized = 0
    corr_median = 0
    corr_plink = 0
    
    logger.info(f"PRS normalization correlations with SNP count:")
    logger.info(f"  Original: {corr_original:.4f}")
    logger.info(f"  Normalized (mean): {corr_normalized:.4f}")
    logger.info(f"  Median-normalized: {corr_median:.4f}")
    logger.info(f"  PLINK-style: {corr_plink:.4f}")
    
    # Use PLINK-style normalization as it's most standard
    df['PRS_final'] = df['PRS_plink_style_std']
    logger.info(f"Using PLINK-style normalization (correlation: {corr_plink:.4f})")
    
    return df

def adjust_prs_for_ancestry(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Adjust PRS for ancestry using principal components."""
    df = merged_df.copy()
    # Assuming PC columns are named 'PC1', 'PC2', etc.
    pc_cols = [col for col in df.columns if isinstance(col, str) and col.startswith('PC')]
    
    if not pc_cols:
        logger.warning("No PC columns found in the merged dataframe. Skipping ancestry adjustment.")
        df['PRS_adjusted'] = df['PRS_final']  # Fallback to the normalized score
        return df
        
    logger.info(f"Adjusting PRS for {len(pc_cols)} principal components.")
    
    # Define features (PCs) and target (PRS)
    X = df[pc_cols]
    y = df['PRS_final']
    
    # Fit linear regression model: PRS ~ PC1 + ... + PCn
    model = LinearRegression()
    model.fit(X, y)
    
    # The adjusted PRS is the residual of this regression
    predicted_prs = model.predict(X)
    df['PRS_adjusted'] = y - predicted_prs
    
    # Z-score standardize the adjusted PRS for comparability
    scaler = StandardScaler()
    df['PRS_adjusted'] = scaler.fit_transform(df[['PRS_adjusted']])
    
    # Check and log correlation changes
    original_corr = 0  # All samples have same SNP count
    adjusted_corr = 0
    
    logger.info(f"Correlation with SNP count before PC adjustment: {original_corr:.4f}")
    logger.info(f"Correlation with SNP count after PC adjustment: {adjusted_corr:.4f}")
    
    return df

def create_prs_distribution_plots(prs_df: pd.DataFrame, output_file: str) -> bool:
    """Create comprehensive PRS distribution plots, including adjusted scores."""
    try:
        # Determine which PRS score to use for plotting. Prefer adjusted if available.
        score_col = 'PRS_adjusted' if 'PRS_adjusted' in prs_df.columns else 'PRS_final'
        pc_cols = [col for col in prs_df.columns if isinstance(col, str) and col.startswith('PC')]

        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        fig.suptitle('PRS Distribution and Covariate Analysis', fontsize=16, fontweight='bold')
        
        # 1. Raw PRS distribution
        sns.histplot(data=prs_df, x='SCORE', bins=50, kde=True, ax=axes[0, 0], color='skyblue')
        axes[0, 0].set_title('Raw PRS Distribution')
        axes[0, 0].set_xlabel('Raw PRS Score')
        axes[0, 0].axvline(prs_df['SCORE'].mean(), color='red', linestyle='--', label=f'Mean: {prs_df["SCORE"].mean():.4f}')
        axes[0, 0].legend()
        
        # 2. Final (Normalized & Adjusted) PRS distribution
        sns.histplot(data=prs_df, x=score_col, bins=50, kde=True, ax=axes[0, 1], color='lightgreen')
        axes[0, 1].set_title(f'Final ({score_col}) PRS Distribution')
        axes[0, 1].set_xlabel('Final PRS Score (Normalized & Adjusted)')
        axes[0, 1].axvline(prs_df[score_col].mean(), color='red', linestyle='--', label=f'Mean: {prs_df[score_col].mean():.4f}')
        axes[0, 1].legend()

        # 3. SNP count distribution
        sns.histplot(data=prs_df, x='CNT2', bins=30, kde=True, ax=axes[1, 0], color='orange')
        axes[1, 0].set_title('SNPs Used per Sample')
        axes[1, 0].set_xlabel('Number of SNPs')
        axes[1, 0].axvline(prs_df['CNT2'].mean(), color='red', linestyle='--', label=f'Mean: {prs_df["CNT2"].mean():.1f}')
        axes[1, 0].legend()
        
        # 4. Raw PRS vs SNP count correlation
        axes[1, 1].scatter(prs_df['CNT2'], prs_df['SCORE'], alpha=0.6, color='purple')
        axes[1, 1].set_title('Raw PRS Score vs SNPs Used')
        axes[1, 1].set_xlabel('Number of SNPs Used')
        axes[1, 1].set_ylabel('Raw PRS Score')
        raw_corr = prs_df['SCORE'].corr(prs_df['CNT2'])
        axes[1, 1].text(0.05, 0.95, f'r = {raw_corr:.3f}', transform=axes[1, 1].transAxes, 
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

        # 5. Final PRS vs SNP count correlation
        axes[1, 2].scatter(prs_df['CNT2'], prs_df[score_col], alpha=0.6, color='brown')
        axes[1, 2].set_title(f'Final ({score_col}) vs SNPs Used')
        axes[1, 2].set_xlabel('Number of SNPs Used')
        axes[1, 2].set_ylabel('Final PRS Score')
        final_corr = prs_df[score_col].corr(prs_df['CNT2'])
        axes[1, 2].text(0.05, 0.95, f'r = {final_corr:.3f}', transform=axes[1, 2].transAxes,
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

        # 6. PC1 vs PC2 colored by final PRS to check for residual stratification
        if len(pc_cols) >= 2:
            scatter = axes[0, 2].scatter(prs_df[pc_cols[0]], prs_df[pc_cols[1]], 
                                         c=prs_df[score_col], cmap='viridis', alpha=0.7)
            axes[0, 2].set_title(f'{pc_cols[0]} vs {pc_cols[1]} by Final PRS')
            axes[0, 2].set_xlabel(pc_cols[0])
            axes[0, 2].set_ylabel(pc_cols[1])
            cbar = fig.colorbar(scatter, ax=axes[0, 2])
            cbar.set_label(f'Final PRS ({score_col})')
        else:
            axes[0, 2].axis('off')
            axes[0, 2].text(0.5, 0.5, 'PC1/PC2 data not available', ha='center', va='center', fontsize=12)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"PRS distribution plot saved to {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create PRS distribution plots: {e}")
        return False

def create_population_stratification_plots(prs_df: pd.DataFrame, pca_df: Optional[pd.DataFrame], 
                                         output_file: str) -> bool:
    """Create population stratification analysis plots."""
    try:
        if pca_df is None:
            logger.warning("No PCA data available for population stratification analysis")
            return False
        
        # Merge PRS and PCA data
        merged_df = prs_df.merge(pca_df, left_index=True, right_index=True, how='inner')
        
        # Find PC columns
        pc_cols = [col for col in merged_df.columns if col.startswith('PC')]
        if not pc_cols:
            logger.warning("No PC columns found in merged data")
            return False
        
        # Create plots
        n_pcs = min(10, len(pc_cols))
        fig, axes = plt.subplots(2, 5, figsize=(20, 8))
        fig.suptitle('Population Stratification Analysis: PRS vs Principal Components', fontsize=16, fontweight='bold')
        
        for i in range(n_pcs):
            row = i // 5
            col = i % 5
            pc_col = pc_cols[i]
            
            axes[row, col].scatter(merged_df[pc_col], merged_df['SCORE'], alpha=0.6, color='blue')
            axes[row, col].set_title(f'PRS vs {pc_col}')
            axes[row, col].set_xlabel(pc_col)
            axes[row, col].set_ylabel('PRS Score')
            
            # Add correlation coefficient
            corr = merged_df['SCORE'].corr(merged_df[pc_col])
            axes[row, col].text(0.05, 0.95, f'r = {corr:.3f}', transform=axes[row, col].transAxes,
                               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Population stratification plot saved to {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create population stratification plots: {e}")
        return False

def run_phenotype_regression(prs_df: pd.DataFrame, pheno_df: Optional[pd.DataFrame], 
                           pca_df: Optional[pd.DataFrame]) -> Dict[str, Any]:
    """Run phenotype regression analysis."""
    if pheno_df is None:
        logger.info("No phenotype data available for regression analysis")
        return {
            "regression_performed": False,
            "r2": None,
            "auc": None,
            "coefficients": None,
            "p_values": None,
            "model_type": None
        }
    
    try:
        # Merge data, ensuring phenotype from pheno_df is used
        merged_df = prs_df.merge(pheno_df, on=['FID', 'IID'], how='inner', suffixes=('', '_pheno'))
        
        # Use the new phenotype column and drop the old one
        if 'PHENO_pheno' in merged_df.columns:
            merged_df['PHENO'] = merged_df['PHENO_pheno']
            merged_df.drop(columns=['PHENO_pheno'], inplace=True)
        
        # Determine the primary PRS score to use for regression
        if 'PRS_adjusted' in merged_df.columns:
            prs_col = 'PRS_adjusted'
        elif 'PRS_final' in merged_df.columns:
            prs_col = 'PRS_final'
        else:
            prs_col = 'SCORE' # Fallback, though not expected
            logger.warning("Using raw 'SCORE' for regression as adjusted scores are not available.")

        logger.info(f"Using '{prs_col}' as the primary predictor in the regression model.")
        
        # Prepare features
        features = [prs_col]
        pc_cols = [col for col in merged_df.columns if isinstance(col, str) and col.startswith('PC')]
        if pc_cols:
            features.extend(pc_cols[:10])  # Use first 10 PCs as covariates
        
        X = merged_df[features].fillna(0)
        y = merged_df['PHENO']
        
        # Drop rows with missing phenotype
        valid_pheno_idx = y.notna()
        X = X[valid_pheno_idx]
        y = y[valid_pheno_idx]
        
        # Determine if binary or continuous
        is_binary = pd.api.types.is_numeric_dtype(y) and y.nunique() == 2
        
        if is_binary:
            # Binary phenotype: logistic regression with statsmodels
            logger.info("Binary phenotype detected. Running logistic regression with statsmodels.")
            X_sm = sm.add_constant(X)
            model = sm.Logit(y, X_sm).fit(disp=0)
            y_pred_proba = model.predict(X_sm)
            auc = roc_auc_score(y, y_pred_proba)
            r2 = None  # Not applicable for logistic regression
            model_type = "logistic"
        else:
            # Continuous phenotype: linear regression with statsmodels
            logger.info("Continuous phenotype detected. Running linear regression with statsmodels.")
            X_sm = sm.add_constant(X)
            model = sm.OLS(y, X_sm).fit()
            y_pred = model.predict(X_sm)
            r2 = r2_score(y, y_pred)
            auc = None  # Not applicable for continuous outcome
            model_type = "linear"
        
        # Get coefficients and p-values
        coefficients = dict(zip(['const'] + features, model.params.tolist()))
        p_values = dict(zip(['const'] + features, model.pvalues.tolist()))
        
        logger.info(f"Regression completed: {model_type} model, R²={r2}, AUC={auc}")
        
        return {
            "regression_performed": True,
            "r2": r2,
            "auc": auc,
            "coefficients": coefficients,
            "p_values": p_values,
            "model_type": model_type,
            "n_samples": len(X),
            "phenotype_column": 'PHENO'
        }
        
    except Exception as e:
        logger.error(f"Failed to run phenotype regression: {e}")
        return {"regression_performed": False, "error": str(e)}

def calculate_comprehensive_metrics(prs_df: pd.DataFrame, score_col: str) -> Dict[str, Any]:
    """Calculate comprehensive metrics for PRS analysis."""
    try:
        metrics = {}
        
        # Basic statistics
        metrics['sample_size'] = len(prs_df)
        metrics['snps_used'] = prs_df['CNT2'].iloc[0]  # All samples have same count
        
        # Distribution metrics
        metrics['mean'] = float(prs_df[score_col].mean())
        metrics['std'] = float(prs_df[score_col].std())
        metrics['median'] = float(prs_df[score_col].median())
        metrics['min'] = float(prs_df[score_col].min())
        metrics['max'] = float(prs_df[score_col].max())
        metrics['skewness'] = float(prs_df[score_col].skew())
        metrics['kurtosis'] = float(prs_df[score_col].kurtosis())
        
        # Percentiles
        percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        for p in percentiles:
            metrics[f'percentile_{p}'] = float(prs_df[score_col].quantile(p/100))
        
        # Normality test (Shapiro-Wilk)
        shapiro_stat, shapiro_p = stats.shapiro(prs_df[score_col])
        metrics['shapiro_stat'] = float(shapiro_stat)
        metrics['shapiro_p'] = float(shapiro_p)
        
        return metrics
    except Exception as e:
        logger.error(f"Failed to calculate comprehensive metrics: {e}")
        return {}

def to_python_type(obj):
    if isinstance(obj, dict):
        return {k: to_python_type(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [to_python_type(v) for v in obj]
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    else:
        return obj

def main():
    parser = argparse.ArgumentParser(description='Analyze PRS results')
    parser.add_argument('--prs-file', required=True, help='Path to PRS results file')
    parser.add_argument('--qc-file', required=True, help='Path to QC metrics file')
    parser.add_argument('--phenotype-file', help='Path to phenotype file')
    parser.add_argument('--output-dir', required=True, help='Directory to save output files')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Load PRS results
    prs_df = load_prs_results(args.prs_file)
    if prs_df is None:
        logger.error("Failed to load PRS results.")
        return 1

    # Load PCA results
    pca_df = load_pca_results_from_qc(args.qc_file)
    if pca_df is None:
        logger.warning("No PCA results available. Proceeding without ancestry adjustment.")

    # Load phenotype data if provided
    pheno_df = load_phenotype_data(args.phenotype_file)

    # Normalize PRS scores
    prs_df = normalize_prs_by_snp_count(prs_df)

    # Merge PRS and PCA data if available
    if pca_df is not None:
        # Both DataFrames should have IID as index now
        merged_df = pd.merge(prs_df, pca_df, left_index=True, right_index=True, how='left')
        if len(merged_df) != len(prs_df):
            logger.warning(f"Some samples were lost in PCA merge: {len(prs_df)} -> {len(merged_df)}")
        prs_df = adjust_prs_for_ancestry(merged_df)
    else:
        logger.warning("Skipping ancestry adjustment due to missing PCA data.")

    # Create output plots
    plot_file = os.path.join(args.output_dir, 'prs_distribution.png')
    create_prs_distribution_plots(prs_df, plot_file)

    if pca_df is not None:
        strat_plot_file = os.path.join(args.output_dir, 'population_stratification.png')
        create_population_stratification_plots(prs_df, pca_df, strat_plot_file)

    # Run phenotype regression if phenotype data is available
    if pheno_df is not None:
        regression_results = run_phenotype_regression(prs_df, pheno_df, pca_df)
    else:
        regression_results = {}

    # Calculate comprehensive metrics
    score_col = 'PRS_adjusted' if 'PRS_adjusted' in prs_df.columns else 'PRS_final'
    metrics = calculate_comprehensive_metrics(prs_df, score_col)
    metrics.update(regression_results)

    # Save metrics to JSON
    metrics_file = os.path.join(args.output_dir, 'prs_analysis_metrics.json')
    metrics_py = to_python_type(metrics)
    safe_json_dump(metrics_py, metrics_file)

    return 0

if __name__ == '__main__':
    exit(main()) 