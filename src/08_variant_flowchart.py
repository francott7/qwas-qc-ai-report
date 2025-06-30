#!/usr/bin/env python3
"""
Generate a flowchart of variant counts after each major filter in the PRS pipeline.
Author: RFtt7

This script automatically collects variant counts from pipeline results and generates
a Mermaid flowchart showing the flow of variants through the pipeline.

Inputs:
    --qc-metrics: QC metrics JSON file (e.g., results/01_qc/qc_metrics.json)
    --effect-sizes: Harmonized effect sizes file (e.g., data/processed/gwas_effect_sizes_filtered_harmonized.txt)
    --causal-snps: Causal SNPs list file (e.g., data/processed/causal_snps.list)
    --output: Output Mermaid flowchart file (e.g., results/08_variant_flowchart.mmd)

Outputs:
    - variant_counts.json: JSON file with step names and variant counts
    - Mermaid flowchart file (e.g., results/08_variant_flowchart.mmd)

Example usage:
    python src/08_variant_flowchart.py \
        --qc-metrics results/01_qc/qc_metrics.json \
        --effect-sizes data/processed/gwas_effect_sizes_filtered_harmonized.txt \
        --causal-snps data/processed/causal_snps.list \
        --output results/08_variant_flowchart.mmd
"""
import argparse
import json
import csv
import sys
import subprocess
from pathlib import Path
from utils import validate_file_exists

def count_lines_in_file(file_path: str) -> int:
    """Count the number of lines in a file."""
    try:
        result = subprocess.run(['wc', '-l', file_path], capture_output=True, text=True)
        if result.returncode == 0:
            return int(result.stdout.strip().split()[0])
        else:
            return 0
    except Exception:
        return 0

def collect_variant_counts(qc_metrics_file: str, effect_sizes_file: str, causal_snps_file: str) -> dict:
    """Collect variant counts from pipeline results."""
    counts = {}
    
    # Raw VCF and After QC counts from QC metrics
    try:
        with open(qc_metrics_file, 'r') as f:
            qc_data = json.load(f)
            counts["Raw VCF"] = qc_data.get("number_of_variants", 0)
            counts["After QC"] = qc_data.get("number_of_variants", 0)  # No filtering in QC step
    except Exception as e:
        print(f"Warning: Could not read QC metrics: {e}")
        counts["Raw VCF"] = 0
        counts["After QC"] = 0
    
    # After Filtering (harmonized effect sizes)
    counts["After Filtering"] = count_lines_in_file(effect_sizes_file)
    
    # After Clumping (same as filtering in this case, no clumping applied)
    counts["After Clumping"] = counts["After Filtering"]
    
    # Causal SNPs
    counts["Causal SNPs"] = count_lines_in_file(causal_snps_file)
    
    # PRS Input (same as causal SNPs)
    counts["PRS Input"] = counts["Causal SNPs"]
    
    return counts

def save_variant_counts(counts: dict, output_file: str):
    """Save variant counts to JSON file."""
    try:
        with open(output_file, 'w') as f:
            json.dump(counts, f, indent=2)
        print(f"Variant counts saved to {output_file}")
    except Exception as e:
        print(f"Error saving variant counts: {e}")

def generate_mermaid(counts_dict):
    """Generate Mermaid flowchart from variant counts."""
    steps = list(counts_dict.keys())
    lines = ["flowchart LR"]
    for i, step in enumerate(steps):
        label = f"{step}\\n{counts_dict[step]:,} variants"
        lines.append(f"  S{i}[\"{label}\"]")
    for i in range(len(steps)-1):
        lines.append(f"  S{i} --> S{i+1}")
    return '\n'.join(lines)

def main():
    parser = argparse.ArgumentParser(description='Generate a variant count flowchart (Mermaid) from pipeline results.')
    parser.add_argument('--qc-metrics', required=True, help='QC metrics JSON file')
    parser.add_argument('--effect-sizes', required=True, help='Harmonized effect sizes file')
    parser.add_argument('--causal-snps', required=True, help='Causal SNPs list file')
    parser.add_argument('--output', required=True, help='Output Mermaid file')
    args = parser.parse_args()

    # Validate input files
    validate_file_exists(args.qc_metrics, "QC metrics")
    validate_file_exists(args.effect_sizes, "Effect sizes")
    validate_file_exists(args.causal_snps, "Causal SNPs")

    # Collect variant counts from pipeline results
    counts = collect_variant_counts(args.qc_metrics, args.effect_sizes, args.causal_snps)
    
    # Save variant counts to JSON
    variant_counts_file = Path(args.output).parent / "variant_counts.json"
    save_variant_counts(counts, str(variant_counts_file))
    
    # Generate and save Mermaid flowchart
    mermaid = generate_mermaid(counts)
    with open(args.output, 'w') as f:
        f.write(mermaid + '\n')
    print(f"Mermaid flowchart saved to {args.output}")

    # Render Mermaid to PNG using Mermaid CLI (mmdc)
    try:
        png_output = str(Path(args.output).parent / "variant_flowchart.png")
        result = subprocess.run([
            "mmdc",
            "-i", args.output,
            "-o", png_output,
            "-b", "transparent",
            "-w", "800",
            "-H", "400"
        ], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Mermaid PNG rendered to {png_output}")
        else:
            print("[Warning] Could not render Mermaid PNG. Is Mermaid CLI (mmdc) installed?")
            print("To install: npm install -g @mermaid-js/mermaid-cli")
            print(result.stderr)
    except FileNotFoundError:
        print("[Warning] Mermaid CLI (mmdc) not found. To install: npm install -g @mermaid-js/mermaid-cli")
    except Exception as e:
        print(f"[Warning] Failed to render Mermaid PNG: {e}")

if __name__ == '__main__':
    main() 