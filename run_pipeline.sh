#!/bin/bash
set -euo pipefail

# GWAS-QC-AI-report: Full pipeline runner
# Usage: bash run_pipeline.sh

# Paths (edit as needed for your data)
VCF=data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
EFFECT_SIZES=data/reference/effect_sizes.csv
CAUSAL_SNPS=data/processed/causal_snps.list
QC_METRICS=results/01_qc/qc_metrics.json
PRS_METRICS=results/04_prs/prs_metrics.json

# 1. QC
echo "[1/9] Running QC..."
python src/01_qc_hail.py --input "$VCF" --output results/01_qc/qc_metrics.json

# 2. QC plots
echo "[2/9] Generating QC plots..."
python src/02_qc_plots.py --qc-metrics results/01_qc/qc_metrics.json --output-dir results/02_qc_plots

# 3. Clumping
echo "[3/9] Running clumping..."
bash src/03_clump_snps.sh

# 4. PRS calculation
echo "[4/9] Calculating PRS..."
python src/04_calculate_prs.py --input "$VCF" --effect-sizes "$EFFECT_SIZES" --output results/04_prs/prs_metrics.json

# 5. PRS analysis
echo "[5/9] Analyzing PRS..."
python src/05_analyze_prs.py --prs-file results/04_prs/prs_metrics.json --qc-file results/01_qc/qc_metrics.json --output-dir results/05_analysis

# 6. Select causal SNPs (optional)
echo "[6/9] Selecting causal SNPs..."
python src/06_select_causal_snps.py --input "$EFFECT_SIZES" --output "$CAUSAL_SNPS"

# 7. Simulate phenotype (optional)
echo "[7/9] Simulating phenotype..."
python src/07_simulate_phenotype.py --causal-snps "$CAUSAL_SNPS" --output data/processed/simulated_phenotype.txt

# 8. Variant flowchart (auto-renders PNG)
echo "[8/9] Generating variant flowchart..."
python src/08_variant_flowchart.py --qc-metrics results/01_qc/qc_metrics.json --effect-sizes data/processed/gwas_effect_sizes_filtered_harmonized.txt --causal-snps "$CAUSAL_SNPS" --output results/08_variant_flowchart/variant_flowchart.mmd

# 9. Generate report (AI-powered)
echo "[9/9] Generating final report..."
python src/09_report.py --qc-metrics results/01_qc/qc_metrics.json --prs-metrics results/04_prs/prs_metrics.json --output-dir reports --template-dir reports_template

echo "\nPipeline complete! See the 'reports/' directory for your final report."
