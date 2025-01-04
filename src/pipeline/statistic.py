import pandas as pd
from cyvcf2 import VCF
import gzip
from collections import defaultdict
from helper.file_utils import save_results_to_csv
from helper.path_define import ground_truth_vcf, statistic_outdir, statistic_nipt_outdir, glimpse_vcf, basevar_vcf, samid
from helper.config import PATHS, PARAMETERS

def process_vcf(vcf_path, method_name):
    """
    Đọc VCF và trích xuất thông tin cần thiết.
    """
    variants = []
    with gzip.open(vcf_path, "rt") as f:
        vcf_reader = VCF.Reader(f)
        for record in vcf_reader:
            variants.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT": str(record.ALT[0]) if record.ALT else ".",
                method_name: True
            })
    return pd.DataFrame(variants)

def compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir):
    """
    So sánh biến thể giữa các phương pháp và lưu kết quả ra file CSV.
    """
    # Load variants from VCF files
    ground_truth_df = process_vcf(ground_truth_path, "GroundTruth")
    basevar_df = process_vcf(basevar_path, "BaseVar")
    glimpse_df = process_vcf(glimpse_path, "Glimpse")

    # Merge all dataframes
    merged_df = pd.merge(ground_truth_df, basevar_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
    merged_df = pd.merge(merged_df, glimpse_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")

    # Fill missing values with False (indicating absence of variant)
    merged_df.fillna(False, inplace=True)

    # Save detailed comparison
    save_results_to_csv("variant_comparison", merged_df, output_dir)

    # Generate summary statistics
    summary = defaultdict(int)
    summary["Total_Variants_GroundTruth"] = ground_truth_df.shape[0]
    summary["Total_Variants_BaseVar"] = basevar_df.shape[0]
    summary["Total_Variants_Glimpse"] = glimpse_df.shape[0]
    summary["Common_Variants"] = merged_df[
        merged_df["GroundTruth"] & merged_df["BaseVar"] & merged_df["Glimpse"]
    ].shape[0]
    summary["Unique_to_GroundTruth"] = merged_df[
        merged_df["GroundTruth"] & ~merged_df["BaseVar"] & ~merged_df["Glimpse"]
    ].shape[0]
    summary["Unique_to_BaseVar"] = merged_df[
        merged_df["BaseVar"] & ~merged_df["GroundTruth"] & ~merged_df["Glimpse"]
    ].shape[0]
    summary["Unique_to_Glimpse"] = merged_df[
        merged_df["Glimpse"] & ~merged_df["GroundTruth"] & ~merged_df["BaseVar"]
    ].shape[0]

    # Save summary statistics
    save_results_to_csv("summary_statistics", pd.DataFrame([summary]), output_dir)

    print(f"Detailed variant comparison and summary statistics saved to {output_dir}")


def run_statistic(fq, chromosome):
    ground_truth_path = ground_truth_vcf(fq, chromosome)
    basevar_path = basevar_vcf(fq, chromosome)
    glimpse_path = glimpse_vcf(fq, chromosome)
    output_dir = statistic_outdir(fq, chromosome)
    compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir)



def statistic(fq, chromosome):
    sample_name = samid(fq)

    if "_" not in sample_name:  # Mẫu đơn
        ground_truth_path = ground_truth_vcf(sample_name, chromosome)
        basevar_path = basevar_vcf(fq, chromosome)
        glimpse_path = glimpse_vcf(fq, chromosome)
        output_dir = statistic_outdir(fq, chromosome)
        compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir)

    else:  # Mẫu NIPT
        child, mom = sample_name.split("_")

        # So sánh với mẹ
        ground_truth_path_mom = ground_truth_vcf(mom, chromosome)
        basevar_path = basevar_vcf(fq, chromosome)
        glimpse_path = glimpse_vcf(fq, chromosome)
        output_dir_mom = statistic_nipt_outdir(fq, chromosome + "mom")
        compare_variants(ground_truth_path_mom, basevar_path, glimpse_path, output_dir_mom)

        # So sánh với con
        ground_truth_path_child = ground_truth_vcf(child, chromosome)
        output_dir_child = statistic_nipt_outdir(fq, chromosome, "child")
        compare_variants(ground_truth_path_child, basevar_path, glimpse_path, output_dir_child)

def run_statistic(fq):
    for chromosome in PARAMETERS["chrs"]:
        statistic(fq, chromosome)