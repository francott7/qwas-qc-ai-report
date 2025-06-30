"""
GWAS QC Report Generation Pipeline
Author: RFtt7

Generate comprehensive GWAS Quality Control and PRS analysis reports using Jinja2 templates
and AI-powered commentary via OpenAI API.

Inputs:
    --qc-metrics: QC metrics JSON file (e.g., results/01_qc/qc_metrics.json)
    --prs-metrics: PRS metrics JSON file (e.g., results/04_prs/prs_metrics.json)
    --template-dir: Directory containing report templates (default: reports_template)
    --output-dir: Output directory for reports (e.g., reports)

Outputs:
    - report_YYYYMMDDHHMM.md: Markdown report with timestamp
    - report_YYYYMMDDHHMM.pdf: PDF report with timestamp (if pandoc is available)

Requirements:
    - config.json with OpenAI API key for AI commentary
    - reports_template/ directory with Jinja2 templates
    - pandoc (optional, for PDF generation)

Example usage:
    python src/09_report.py \
        --qc-metrics results/01_qc/qc_metrics.json \
        --prs-metrics results/04_prs/prs_metrics.json \
        --template-dir reports_template \
        --output-dir reports
"""

#!/usr/bin/env python3

import argparse
import json
import os
from typing import Dict, Any
from jinja2 import Environment, FileSystemLoader
from openai import OpenAI
import subprocess
from datetime import datetime
from pathlib import Path
import yaml

# Load API key from config file
config_path = Path('config.json')
if config_path.exists():
    with open(config_path) as f:
        config = json.load(f)
        os.environ['OPENAI_API_KEY'] = config.get('openai_api_key', '')

def load_metrics(metrics_files: Dict[str, str]) -> Dict[str, Any]:
    """Load all metrics from JSON files."""
    metrics = {}
    for key, file_path in metrics_files.items():
        with open(file_path, 'r') as f:
            metrics[key] = json.load(f)
    return metrics

def load_prs_analysis_metrics(path: str) -> dict:
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except Exception:
        return {}

def load_clumping_summary(path: str) -> dict:
    summary = {
        'clump_p1': None,
        'clump_r2': None,
        'clump_kb': None,
        'variants_before': None,
        'variants_after': None,
        'warnings': [],
        'notes': []
    }
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if '--clump-p1' in line:
                summary['clump_p1'] = line.split('--clump-p1')[1].strip().split()[0]
            if '--clump-r2' in line:
                summary['clump_r2'] = line.split('--clump-r2')[1].strip().split()[0]
            if '--clump-kb' in line:
                summary['clump_kb'] = line.split('--clump-kb')[1].strip().split()[0]
            if 'variants loaded from .bim file' in line:
                summary['variants_before'] = int(line.split('variants loaded')[0].strip().split()[-1])
            if 'variants and' in line and 'pass filters and QC' in line:
                summary['variants_after'] = int(line.split('variants and')[0].strip().split()[-1])
            if 'Warning:' in line:
                summary['warnings'].append(line.strip())
            if 'Note:' in line:
                summary['notes'].append(line.strip())
        return summary
    except Exception:
        return summary

def load_flowchart(path: str) -> str:
    try:
        with open(path, 'r') as f:
            return f.read()
    except Exception:
        return "Flowchart not found."

def generate_ai_commentary(metrics: Dict[str, Any]) -> str:
    """Generate AI commentary for the report using OpenAI API."""
    # Load API key from config file
    with open('config.json', 'r') as f:
        config = json.load(f)
        api_key = config.get('openai_api_key')
    
    if not api_key:
        return "API key not found in config.json"
    
    # Initialize OpenAI client
    client = OpenAI(api_key=api_key)
    
    # Read the prompt template from file
    prompt_path = Path('prompts/report_prompt.txt')
    if not prompt_path.exists():
        return "Prompt template file not found"
    
    with open(prompt_path, 'r') as f:
        prompt_template = f.read()
    
    # Calculate heterozygosity rate safely
    n_het = metrics['qc_metrics']['sample_qc'].get('n_het', 0)
    n_non_ref = metrics['qc_metrics']['sample_qc'].get('n_non_ref', 0)
    heterozygosity_rate = f"{n_het / n_non_ref:.3f}" if n_non_ref > 0 else 'N/A'
    
    # Handle missing gwas_qc section
    gwas_qc = metrics['qc_metrics'].get('gwas_qc', {})
    n_snps = gwas_qc.get('n_snps', 'N/A')
    effect_size_min = f"{gwas_qc.get('effect_size_min', 'N/A'):.3f}" if isinstance(gwas_qc.get('effect_size_min'), (int, float)) else 'N/A'
    effect_size_max = f"{gwas_qc.get('effect_size_max', 'N/A'):.3f}" if isinstance(gwas_qc.get('effect_size_max'), (int, float)) else 'N/A'
    effect_size_mean = f"{gwas_qc.get('effect_size_mean', 'N/A'):.3f}" if isinstance(gwas_qc.get('effect_size_mean'), (int, float)) else 'N/A'
    
    # Only pass summary statistics for PRS
    prs_metrics = metrics.get('prs_metrics', {})
    formatted_metrics = {
        'sample_call_rate': f"{metrics['qc_metrics']['sample_qc']['call_rate']:.3f}",
        'heterozygosity_rate': heterozygosity_rate,
        'variant_call_rate': f"{metrics['qc_metrics']['variant_qc']['call_rate']:.3f}",
        'maf': f"{metrics['qc_metrics']['variant_qc']['allele_frequency']:.3f}",
        'avg_quality': f"{metrics['qc_metrics']['dosage_qc']['avg_quality']:.3f}",
        'snps_below_threshold': metrics['qc_metrics']['dosage_qc']['snps_below_threshold'],
        'n_snps': n_snps,
        'effect_size_min': effect_size_min,
        'effect_size_max': effect_size_max,
        'effect_size_mean': effect_size_mean,
        'prs_mean': prs_metrics.get('mean', 'N/A'),
        'prs_std': prs_metrics.get('std', 'N/A'),
        'prs_auc': prs_metrics.get('auc', 'N/A'),
        'prs_r2': prs_metrics.get('r2', 'N/A')
    }
    
    # Format the prompt with metrics
    prompt = prompt_template.format(**formatted_metrics)
    
    try:
        # Call OpenAI API
        response = client.chat.completions.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": "You are an expert in GWAS and PRS analysis."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7,
            max_tokens=1000
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error generating commentary: {str(e)}"

def to_yaml(value):
    """Convert a Python object to YAML string."""
    return yaml.dump(value, default_flow_style=False)

def render_report(metrics: Dict[str, Any], template_dir: str, output_dir: str):
    """Render the QC report using templates."""
    # Convert template_dir to absolute path
    template_dir = os.path.abspath(template_dir)
    
    # Set up Jinja2 environment
    env = Environment(loader=FileSystemLoader(template_dir))
    env.filters['to_yaml'] = to_yaml
    template = env.get_template('report.md.j2')
    
    # Generate AI commentary
    commentary = generate_ai_commentary(metrics)
    
    # Add pipeline parameters
    metrics['pipeline_parameters'] = {
        'maf_threshold': 0.01,
        'call_rate_threshold': 0.95,
        'hwe_threshold': 1e-6,
        'dosage_threshold': 0.1
    }

    # Load clumping summary
    clumping_summary = load_clumping_summary('results/03_clumping/chr22_clumped.log')
    # If no variants after clumping, but PRS analysis used a subset, update the count
    if clumping_summary['variants_after'] == clumping_summary['variants_before'] and 'prs_analysis' in metrics and metrics['prs_analysis'].get('snps_used'):
        clumping_summary['variants_after'] = metrics['prs_analysis']['snps_used']

    # Load variant flowchart
    flowchart_code = load_flowchart('results/08_variant_flowchart.mmd')

    # Ensure dosage QC is passed as a dict
    if 'dosage_qc' not in metrics['qc_metrics']:
        metrics['qc_metrics']['dosage_qc'] = {'avg_quality': 0, 'snps_below_threshold': 0, 'distribution': {'0': 0, '1': 0, '2': 0}}
    # Ensure PCA explained variance is a list of floats
    if 'pca' in metrics['qc_metrics'] and 'explained_variance' in metrics['qc_metrics']['pca']:
        metrics['qc_metrics']['pca']['explained_variance'] = [float(x) for x in metrics['qc_metrics']['pca']['explained_variance']]

    # Build metrics dict for template
    metrics = {
        'qc_metrics': metrics['qc_metrics'],
        'prs_metrics': metrics['prs_metrics'],
        'prs_analysis': metrics['prs_analysis'],
        'clumping_summary': clumping_summary,
        'variant_flowchart': flowchart_code,
        'figures': {
            'quality_metrics': 'results/02_qc_plots/quality_metrics.png',
            'dosage_distribution': 'results/02_qc_plots/dosage_distribution.png',
            'pca_plot': 'results/02_qc_plots/pca_plot.png',
            'prs_distribution': 'results/05_analysis/prs_distribution.png',
            'population_stratification': 'results/05_analysis/population_stratification.png',
        },
        'pipeline_parameters': metrics['pipeline_parameters']
    }

    # Check if flowchart PNG exists
    flowchart_img_path = 'results/08_variant_flowchart/variant_flowchart.png'
    flowchart_img_exists = os.path.isfile(flowchart_img_path)

    # Render report
    report = template.render(
        metrics=metrics,
        commentary=commentary,
        date=datetime.now().strftime('%Y-%m-%d'),
        version='1.0',
        flowchart_img_exists=flowchart_img_exists,
        flowchart_img_path=flowchart_img_path
    )
    
    # Create timestamp for unique filename
    timestamp = datetime.now().strftime('%Y%m%d%H%M')
    base_filename = f'report_{timestamp}'
    md_output = os.path.join(output_dir, f'{base_filename}.md')
    pdf_output = os.path.join(output_dir, f'{base_filename}.pdf')
    
    # Save report
    with open(md_output, 'w') as f:
        f.write(report)
    
    # Convert to PDF
    os.system(f'pandoc {md_output} -o {pdf_output} --pdf-engine=xelatex')

def main():
    parser = argparse.ArgumentParser(description='Generate GWAS QC report')
    parser.add_argument('--qc-metrics', required=True,
                       help='QC metrics JSON file')
    parser.add_argument('--prs-metrics', required=True,
                       help='PRS metrics JSON file')
    parser.add_argument('--template-dir', default='reports_template',
                        help='Directory containing report templates')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for reports')
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load metrics
    metrics = load_metrics({
        'qc_metrics': args.qc_metrics,
        'prs_metrics': args.prs_metrics
    })

    # Load PRS analysis metrics
    prs_analysis = load_prs_analysis_metrics('results/05_analysis/prs_analysis_metrics.json')
    metrics['prs_analysis'] = prs_analysis

    # Generate markdown report
    render_report(metrics, args.template_dir, args.output_dir)

if __name__ == '__main__':
    main() 