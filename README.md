# Mito_Plot
Variant Plot Considering All Heterozygosity Rates in Large-Scale Samples

# Mitochondrial Heteroplasmy Analysis and Visualization

This repository contains scripts to calculate mitochondrial DNA heteroplasmy from VCF files and visualize the results in an interactive circular plot.

In addition to the original 2D interactive plots which allow free zooming and panning, the tool now supports generating interactive 3D circular plots.

This enhancement enables detailed visualization of aggregated variants and heteroplasmy rates across large cohorts, providing better insight into complex mutation patterns that were difficult to explore with previous 2D-only visualizations.

## Files

- `step1.py`:  
  Processes multiple mitochondrial VCF files (*.vcf) to calculate heteroplasmy rates at each variant position. Outputs a merged TSV file summarizing heteroplasmy rates per sample and position.

- `step2.py`:  
  Loads the merged TSV from run_step1.py and generates an interactive Plotly circular plot showing heteroplasmy rates by position with mitochondrial gene regions annotated.
  This script creates a 2D figure, where zoom in/out and panning on desired regions are freely possible for detailed inspection.

- `step3.py`:  
  Loads the merged TSV from run_step1.py and generates an interactive 3D circular plot representing heteroplasmy rates across samples.
  The plot positions samples along the Z-axis with compressible height for better visualization of multiple samples and mitochondrial gene regions annotated in 3D.

## Usage

1. Prepare a directory containing your mitochondrial VCF files named like `sample1_chrM.vcf`.

2. Run step Plotting a 2D Plot:
python step2.py All_YYYMMDD.tsv

3. Run step Plotting a 3D Plot:
python step3.py All_YYYMMDD.tsv


[MT_Plot_Sample_2D.pdf](https://github.com/user-attachments/files/22736247/MT_Plot_Sample_2D.pdf)

[MT_Plot_Sample_3D.pdf](https://github.com/user-attachments/files/22736176/MT_Plot_Sample_3D.pdf)

[MT_Plot_Sample_test.html](https://github.com/IDRC-MITO/Mito_Plot/blob/main/Mito_Heteroplasmy_Radial3D_Height_Lines.html)


## Input File Specification

The software expects input VCF files generated on the GRCh38 reference genome (VCF v4 format), with the chromosome column (CHROM) labeled as chrM. This is the standard convention for mitochondrial DNA in GRCh38-based VCF data.
