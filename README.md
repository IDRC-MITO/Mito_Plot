# Mito_Plot
Variant Plot Considering All Heterozygosity Rates in Large-Scale Samples

# Mitochondrial Heteroplasmy Analysis and Visualization

This repository contains scripts to calculate mitochondrial DNA heteroplasmy from VCF files and visualize the results in an interactive circular plot.

## Files

- `run_step1.py`:  
  Processes multiple mitochondrial VCF files (*.vcf) to calculate heteroplasmy rates at each variant position. Outputs a merged TSV file summarizing heteroplasmy rates per sample and position.

- `run_step2.py`:  
  Loads the merged TSV from `run_step1.py` and generates an interactive Plotly circular plot showing heteroplasmy rates by position with mitochondrial gene regions annotated.

## Usage

1. Prepare a directory containing your mitochondrial VCF files named like `sample1_chrM.vcf`.

2. Run step 1 to generate the merged TSV file:
python step2.py All_YYYMMDD.tsv


[MT_Plot.pdf](https://github.com/user-attachments/files/22651014/MT_Plot.pdf)
