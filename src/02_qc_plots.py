#!/usr/bin/env python3

"""
Generate QC plots from GWAS QC metrics.
Author: RFtt7

Inputs:
    --metrics-file: Path to QC metrics JSON file (e.g., results/01_qc/qc_metrics.json)

Outputs:
    - QQ plot (e.g., results/02_qc_plots/qq_plot.png)
    - Manhattan plot (e.g., results/02_qc_plots/manhattan_plot.png)
    - PCA explained variance plot (e.g., results/02_qc_plots/pca_plot.png)
    - Dosage distribution plot (e.g., results/02_qc_plots/dosage_distribution.png)
    - Quality metrics plot (e.g., results/02_qc_plots/quality_metrics.png)

Example usage:
    python src/02_qc_plots.py \
        --metrics-file results/01_qc/qc_metrics.json \
        --output-dir results/02_qc_plots
"""

import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import os
from pathlib import Path
import logging
from typing import Dict, Any, List, Tuple, Optional
from utils import validate_file_exists, safe_json_dump

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Input: results/01_qc/qc_metrics.json

def load_metrics(metrics_file: str) -> Optional[Dict[str, Any]]:
    """Load QC metrics from JSON file."""
    if not validate_file_exists(metrics_file, "Metrics"):
        return None
        
    try:
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)
        logger.info("Successfully loaded metrics file")
        return metrics
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing metrics file: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error loading metrics file: {e}")
        return None

def annotate_no_data(ax, message="No Data"):
    ax.text(0.5, 0.5, message, fontsize=18, color='red', ha='center', va='center', alpha=0.7, transform=ax.transAxes)
    ax.set_axis_off()

def create_qq_plot(p_values: List[float], output_file: str) -> bool:
    """Create QQ plot for p-values."""
    try:
        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        if not p_values or all([p == 0 or p is None for p in p_values]):
            annotate_no_data(ax)
            plt.savefig(output_file)
            plt.close()
            logger.warning(f"QQ plot not generated: no p-values provided. Annotated as 'No Data'.")
            return False
        
        # Sort p-values
        p_values = np.sort(p_values)
        
        # Calculate expected p-values
        n = len(p_values)
        expected = -np.log10(np.linspace(1/n, 1, n))
        observed = -np.log10(p_values)
        
        # Plot
        plt.scatter(expected, observed, alpha=0.5)
        plt.plot([0, max(expected)], [0, max(expected)], 'r--')
        
        # Labels and title
        plt.xlabel('Expected -log10(p)')
        plt.ylabel('Observed -log10(p)')
        plt.title('Q-Q Plot')
        
        # Save
        plt.savefig(output_file)
        plt.close()
        logger.info(f"QQ plot saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating QQ plot: {e}")
        return False

def create_manhattan_plot(chromosomes: List[str], positions: List[int], 
                         p_values: List[float], output_file: str) -> bool:
    """Create Manhattan plot."""
    try:
        plt.figure(figsize=(15, 5))
        ax = plt.gca()
        if not chromosomes or not positions or not p_values or all([p == 0 or p is None for p in p_values]):
            annotate_no_data(ax)
            plt.savefig(output_file)
            plt.close()
            logger.warning(f"Manhattan plot not generated: missing or empty data. Annotated as 'No Data'.")
            return False
        
        # Convert chromosomes to numeric values for plotting
        chrom_dict = {chrom: i for i, chrom in enumerate(sorted(set(chromosomes)))}
        chrom_nums = [chrom_dict[chrom] for chrom in chromosomes]
        
        # Plot
        plt.scatter(positions, -np.log10(p_values), c=chrom_nums, alpha=0.5)
        
        # Add chromosome labels
        plt.xticks(range(len(chrom_dict)), sorted(set(chromosomes)), rotation=45)
        
        # Labels and title
        plt.xlabel('Chromosome')
        plt.ylabel('-log10(p)')
        plt.title('Manhattan Plot')
        
        # Save
        plt.savefig(output_file)
        plt.close()
        logger.info(f"Manhattan plot saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating Manhattan plot: {e}")
        return False

def create_pca_plot(pca_metrics: Dict[str, Any], output_file: str) -> bool:
    """Create PCA plot showing explained variance."""
    try:
        plt.figure(figsize=(10, 5))
        ax = plt.gca()
        if not pca_metrics or 'explained_variance' not in pca_metrics or not pca_metrics['explained_variance']:
            annotate_no_data(ax)
            plt.savefig(output_file)
            plt.close()
            logger.warning(f"PCA plot not generated: missing explained variance. Annotated as 'No Data'.")
            return False
        
        # Plot explained variance
        plt.bar(range(len(pca_metrics['explained_variance'])), 
                pca_metrics['explained_variance'])
        
        # Labels and title
        plt.xlabel('Principal Component')
        plt.ylabel('Explained Variance')
        plt.title('PCA Explained Variance')
        
        # Save
        plt.savefig(output_file)
        plt.close()
        logger.info(f"PCA plot saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating PCA plot: {e}")
        return False

def create_dosage_distribution_plot(dosage_metrics: Dict[str, Any], output_file: str) -> bool:
    """Create plot showing dosage distribution."""
    try:
        plt.figure(figsize=(10, 5))
        ax = plt.gca()
        distribution = dosage_metrics.get('distribution', {})
        if not distribution or all([v == 0 for v in distribution.values()]):
            annotate_no_data(ax)
            plt.savefig(output_file)
            plt.close()
            logger.warning(f"Dosage distribution plot not generated: no dosage data. Annotated as 'No Data'.")
            return False
        
        # Extract distribution data
        categories = list(distribution.keys())
        values = list(distribution.values())
        
        # Create bar plot
        plt.bar(categories, values)
        
        # Labels and title
        plt.xlabel('Dosage Category')
        plt.ylabel('Count')
        plt.title('Dosage Distribution')
        
        # Save
        plt.savefig(output_file)
        plt.close()
        logger.info(f"Dosage distribution plot saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating dosage distribution plot: {e}")
        return False

def create_quality_metrics_plot(metrics: Dict[str, Any], output_file: str) -> bool:
    """Create plot showing various quality metrics."""
    try:
        plt.figure(figsize=(12, 6))
        ax = plt.gca()
        metrics_to_plot = {
            'Call Rate': metrics.get('variant_qc', {}).get('call_rate', 0),
            'MAF': metrics.get('variant_qc', {}).get('allele_frequency', 0),
            'Ts/Tv Ratio': metrics.get('ts_tv_ratio', 0),
            'Dosage Quality': metrics.get('dosage_qc', {}).get('avg_quality', 0)
        }
        if all([v == 0 or v is None for v in metrics_to_plot.values()]):
            annotate_no_data(ax)
            plt.tight_layout()
            plt.savefig(output_file)
            plt.close()
            logger.warning(f"Quality metrics plot not generated: all metrics missing or zero. Annotated as 'No Data'.")
            return False
        
        # Create bar plot
        plt.bar(metrics_to_plot.keys(), metrics_to_plot.values())
        
        # Labels and title
        plt.xlabel('Metric')
        plt.ylabel('Value')
        plt.title('Quality Metrics Summary')
        plt.xticks(rotation=45)
        
        # Save
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        logger.info(f"Quality metrics plot saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating quality metrics plot: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Generate QC plots')
    parser.add_argument('--metrics', required=True, help='Input JSON file with QC metrics')
    parser.add_argument('--output-dir', required=False, default='results/plots', help='Output directory for plots (default: results/plots)')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Load metrics
    metrics = load_metrics(args.metrics)
    if not metrics:
        logger.error("Failed to load metrics file")
        exit(1)

    # Generate plots
    success = True
    plot_status = {}
    
    # Create quality metrics plot
    plot_status['quality_metrics'] = create_quality_metrics_plot(metrics, os.path.join(args.output_dir, 'quality_metrics.png'))
    
    # Create dosage distribution plot if dosage metrics are available
    if 'dosage_qc' in metrics:
        plot_status['dosage_distribution'] = create_dosage_distribution_plot(metrics['dosage_qc'], os.path.join(args.output_dir, 'dosage_distribution.png'))
    
    # Create PCA plot if PCA metrics are available
    if 'pca' in metrics:
        plot_status['pca'] = create_pca_plot(metrics['pca'], os.path.join(args.output_dir, 'pca_plot.png'))
    
    # Create QQ and Manhattan plots if p-values are available
    if 'p_values' in metrics:
        plot_status['qq'] = create_qq_plot(metrics['p_values'], os.path.join(args.output_dir, 'qq_plot.png'))
        
        if not create_manhattan_plot(metrics.get('chromosomes', []),
                                   metrics.get('positions', []),
                                   metrics['p_values'],
                                   os.path.join(args.output_dir, 'manhattan_plot.png')):
            success = False
    
    # Print/log summary
    for plot, status in plot_status.items():
        if status:
            logger.info(f"{plot} plot generated successfully.")
        else:
            logger.warning(f"{plot} plot was not generated or is empty.")

    if not success:
        logger.error("Some plots failed to generate")
        exit(1)
    
    logger.info("All plots generated successfully")

if __name__ == '__main__':
    main() 