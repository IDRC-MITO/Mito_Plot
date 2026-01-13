#!/usr/bin/env python3
"""
Mitochondrial Heteroplasmy Calculator from VCF files.

Merges heteroplasmy rates across multiple chrM VCF files into a TSV matrix.
Supports both uncompressed and gzipped VCF files.

Usage:
    python heteroplasmy_calculator.py /path/to/vcf_directory [--output output.tsv]

Arguments:
    input_dir    Directory containing *chrM.vcf or *chrM.vcf.gz files
    --output, -o Output TSV file (default: All_YYYYMMDD.tsv)
"""

import os
import glob
import argparse
import gzip
from datetime import date

def open_vcf(vcf_file):
    """Open VCF file, handling both plain and gzipped formats."""
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    return open(vcf_file, 'r')

def calculate_heteroplasmy_rate(vcf_file):
    """Calculate heteroplasmy rates for variants in a single VCF file."""
    results = {}
    with open_vcf(vcf_file) as file:
        for line in file:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            format_fields = fields[8].split(":")
            sample_data = fields[9].split(":")
            
            if "GT" not in format_fields or "AD" not in format_fields:
                continue
            
            try:
                gt_index = format_fields.index("GT")
                ad_index = format_fields.index("AD")
                
                if len(sample_data) <= max(gt_index, ad_index):
                    continue
                
                gt = sample_data[gt_index]
                ad = sample_data[ad_index]
                
                if "/" in gt or "|" in gt:
                    ad_values = list(map(int, ad.split(",")))
                    total_depth = sum(ad_values)
                    if total_depth > 0:
                        alt_allele_depth = sum(ad_values[1:])
                        heteroplasmy_rate = (alt_allele_depth / total_depth) * 100
                        results[pos] = {"ref": ref, "alt": alt, "rate": f"{heteroplasmy_rate:.2f}%"}
            except (ValueError, ZeroDivisionError):
                continue
    
    return results

def process_all_vcf_files(directory):
    """Process all *chrM.vcf* files in the directory."""
    all_results = {}
    positions = set()

    vcf_files = glob.glob(os.path.join(directory, "*chrM.vcf*"))
    if not vcf_files:
        raise ValueError(f"No *chrM.vcf* files found in {directory}")

    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace("_chrM.vcf", "").replace(".vcf", "").replace(".gz", "")
        print(f"Processing: {vcf_file}")
        try:
            results = calculate_heteroplasmy_rate(vcf_file)
            if not results:
                print(f"Warning: No variants found in {vcf_file}")
            else:
                print(f"Found {len(results)} variants")
            all_results[sample_name] = results
            positions.update(results.keys())
        except Exception as e:
            print(f"Error processing {vcf_file}: {str(e)}")
            continue

    return all_results, sorted(positions)

def save_merged_results(all_results, positions, output_file):
    """Save results as TSV with samples as columns."""
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    with open(output_file, 'w', newline='') as tsvfile:
        header = ["chrM", "Position", "Ref", "Alt"] + sorted(all_results.keys())
        tsvfile.write("\t".join(header) + "\n")
        
        for pos in positions:
            ref = ""
            alt = ""
            row = ["chrM", pos]
            for sample in sorted(all_results.keys()):
                if pos in all_results[sample]:
                    if not ref:
                        ref = all_results[sample][pos]["ref"]
                    if not alt:
                        alt = all_results[sample][pos]["alt"]
                    row.append(all_results[sample][pos]["rate"])
                else:
                    row.append("")
            if ref:
                row.insert(2, ref)
            if alt:
                row.insert(3, alt)
            tsvfile.write("\t".join(row) + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Calculate and merge mitochondrial heteroplasmy rates from VCF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python heteroplasmy_calculator.py ./vcf_files
  python heteroplasmy_calculator.py /path/to/vcf -o results.tsv
        """
    )
    parser.add_argument("input_dir", help="Directory containing *chrM.vcf* files")
    parser.add_argument("-o", "--output", help="Output TSV file (default: All_YYYYMMDD.tsv)")
    args = parser.parse_args()

    output_filename = args.output or f"All_{date.today().strftime('%Y%m%d')}.tsv"
    
    all_results, positions = process_all_vcf_files(args.input_dir)
    save_merged_results(all_results, positions, output_filename)

    print(f"Results saved to '{output_filename}'")
    print(f"Samples processed: {len(all_results)}")
    print(f"Total variants: {len(positions)}")

if __name__ == "__main__":
    main()

