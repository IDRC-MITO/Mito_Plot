import os
import glob
from datetime import date

def calculate_heteroplasmy_rate(vcf_file):
    results = {}
    with open(vcf_file, 'r') as file:
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
            
            gt_index = format_fields.index("GT")
            ad_index = format_fields.index("AD")
            
            if len(sample_data) <= max(gt_index, ad_index):
                continue
            
            gt = sample_data[gt_index]
            ad = sample_data[ad_index]
            
            if "/" in gt or "|" in gt:
                try:
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
    all_results = {}
    positions = set()

    for vcf_file in glob.glob(os.path.join(directory, "*chrM.vcf")):
        sample_name = os.path.basename(vcf_file).replace("_chrM.vcf", "")
        print(f"Processing file: {vcf_file}")
        try:
            results = calculate_heteroplasmy_rate(vcf_file)
            if not results:
                print(f"Warning: No variants found in {vcf_file}")
            else:
                print(f"Found {len(results)} variants in {vcf_file}")
            all_results[sample_name] = results
            positions.update(results.keys())
        except Exception as e:
            print(f"Error processing file {vcf_file}: {str(e)}")
            continue

    return all_results, sorted(positions)

def save_merged_results(all_results, positions, output_file):
    with open(output_file, 'w', newline='') as tsvfile:
        header = ["chrM", "Position", "Ref", "Alt"] + list(all_results.keys())
        tsvfile.write("\t".join(header) + "\n")
        
        for pos in positions:
            ref = ""
            alt = ""
            row = ["chrM", pos]
            for sample in all_results.keys():
                if pos in all_results[sample]:
                    if not ref and not alt:
                        ref = all_results[sample][pos]["ref"]
                        alt = all_results[sample][pos]["alt"]
                    row.append(all_results[sample][pos]["rate"])
                else:
                    row.append("")
            row.insert(2, ref)
            row.insert(3, alt)
            tsvfile.write("\t".join(row) + "\n")


if __name__ == "__main__":
    # Please adjust this directory path to the folder containing your VCF files
    vcf_directory = "./vcf_files"
    output_filename = f"All_{date.today().strftime('%Y%m%d')}.tsv"

    all_results, positions = process_all_vcf_files(vcf_directory)
    save_merged_results(all_results, positions, output_filename)

    print(f"Results saved to '{output_filename}'.")
    print(f"Number of processed samples: {len(all_results)}")
    print(f"Total number of variants: {len(positions)}")
