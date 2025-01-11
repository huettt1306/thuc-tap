import pandas as pd
import os
from cyvcf2 import VCF
from collections import defaultdict
from helper.file_utils import save_results_to_csv
from helper.path_define import ground_truth_vcf, statistic_outdir, statistic_nipt_outdir, glimpse_vcf, basevar_vcf, samid
from helper.config import PATHS, PARAMETERS
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "statistic_pipeline.log"))


def process_vcf(vcf_path, method_name):
    """
    Đọc VCF và trích xuất thông tin cần thiết.
    """
    variants = []
    try:
        logger.info(f"Processing VCF file: {vcf_path} for method {method_name}")
        vcf_reader = VCF(vcf_path)
        for record in vcf_reader:
            variants.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT": str(record.ALT[0]) if record.ALT else ".",
                f"INFO_{method_name}": dict(record.INFO),
                method_name: True
            })
    except Exception as e:
        logger.error(f"Error processing VCF file {vcf_path}: {e}")
        raise

    logger.info(f"Finished processing VCF file: {vcf_path}")
    return pd.DataFrame(variants)

def compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir):
    """
    So sánh biến thể giữa các phương pháp và lưu kết quả ra file CSV.
    """
    logger.info(f"Comparing variants for paths: {ground_truth_path}, {basevar_path}, {glimpse_path}")

    try:
        ground_truth_df = process_vcf(ground_truth_path, "Truth")
        basevar_df = process_vcf(basevar_path, "BaseVar")
        glimpse_df = process_vcf(glimpse_path, "Glimpse")

        merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df = pd.merge(ground_truth_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df.fillna(False, inplace=True)

        save_results_to_csv("variants", merged_df, output_dir)

        logger.info(f"Detailed variant comparison and summary statistics saved to {output_dir}")
    except Exception as e:
        logger.error(f"Error comparing variants: {e}")
        raise

def statistic(fq, chromosome):
    sample_name = samid(fq)

    try:
        logger.info(f"Starting statistics for sample {sample_name} on chromosome {chromosome}")

        if "_" not in sample_name:
            ground_truth_path = ground_truth_vcf(sample_name, chromosome)
            basevar_path = basevar_vcf(fq, chromosome)
            glimpse_path = glimpse_vcf(fq, chromosome)
            output_dir = statistic_outdir(fq, chromosome)

            compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir)
        else:
            child, mom = sample_name.split("_")

            ground_truth_path_mom = ground_truth_vcf(mom, chromosome)
            output_dir_mom = statistic_nipt_outdir(fq, chromosome, "mom")
            compare_variants(ground_truth_path_mom, basevar_vcf(fq, chromosome), glimpse_vcf(fq, chromosome), output_dir_mom)

            ground_truth_path_child = ground_truth_vcf(child, chromosome)
            output_dir_child = statistic_nipt_outdir(fq, chromosome, "child")
            compare_variants(ground_truth_path_child, basevar_vcf(fq, chromosome), glimpse_vcf(fq, chromosome), output_dir_child)

        logger.info(f"Completed statistics for sample {sample_name} on chromosome {chromosome}")
    except Exception as e:
        logger.error(f"Error in statistic function for sample {sample_name} and chromosome {chromosome}: {e}")
        raise

def run_statistic(fq):
    try:
        logger.info(f"Starting full statistics pipeline for sample {fq}")

        for chromosome in PARAMETERS["chrs"]:
            statistic(fq, chromosome)

        logger.info(f"Completed full statistics pipeline for sample {fq}")
    except Exception as e:
        logger.error(f"Error in run_statistic pipeline for sample {fq}: {e}")
        raise
