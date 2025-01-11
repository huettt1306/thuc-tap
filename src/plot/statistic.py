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

def is_rare_variant(info, key="AF_EAS", threshold=0.01):
    """
    Kiểm tra xem biến thể có phải là biến thể hiếm (AF_ESA < threshold) hay không.
    """
    try:
        af_esa = float(next((entry.split("=")[1] for entry in info.split(";") if entry.startswith(key)), "1"))
        return af_esa < threshold
    except ValueError:
        return False



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
                "INFO": str(record.INFO) if method_name == "Truth" else None,
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

        merged_df["Truth"] = merged_df["Truth"].fillna(False)
        merged_df["BaseVar"] = merged_df["BaseVar"].fillna(False)
        merged_df["Glimpse"] = merged_df["Glimpse"].fillna(False)
        merged_df["INFO"] = merged_df["INFO"].fillna("")

        save_results_to_csv("variants", merged_df, output_dir)
        logger.info(f"Detailed variant comparison and summary statistics saved to {output_dir}")

        return merged_df
    except Exception as e:
        logger.error(f"Error comparing variants: {e}")
        raise


def calculate_statistics(df, label="Total"):
    """
    Tính các chỉ số thống kê cho một DataFrame và trả về kết quả dưới dạng DataFrame.
    """
    def calculate_metrics(truth_col, method_col):
        true_positive = len(df[df[truth_col] & df[method_col]])
        false_positive = len(df[~df[truth_col] & df[method_col]])
        false_negative = len(df[df[truth_col] & ~df[method_col]])
        true_negative = len(df[~df[truth_col] & ~df[method_col]])

        accuracy = (true_positive + true_negative) / len(df) if len(df) > 0 else 0
        recall = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0
        precision = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0
        f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        return {
            "True Positives": true_positive,
            "False Positives": false_positive,
            "False Negatives": false_negative,
            "True Negatives": true_negative,
            "Accuracy": accuracy,
            "Recall": recall,
            "F1 Score": f1_score,
        }

    # Tổng số biến thể
    total_variants = len(df)
    truth_variants = len(df[df["Truth"]])
    basevar_variants = len(df[df["BaseVar"]])
    glimpse_variants = len(df[df["Glimpse"]])

    # Tính chỉ số cho BaseVar và Glimpse
    basevar_metrics = calculate_metrics("Truth", "BaseVar")
    glimpse_metrics = calculate_metrics("Truth", "Glimpse")

    # Tạo DataFrame với các chỉ số thống kê
    stats_data = {
        "Category": [label] * 2,
        "Method": ["BaseVar", "Glimpse"],
        "Total Variants": [total_variants, total_variants],
        "Truth Variants": [truth_variants, truth_variants],
        "BaseVar Variants": [basevar_variants, ""],  # Chỉ cho BaseVar
        "Glimpse Variants": ["", glimpse_variants],  # Chỉ cho Glimpse
        "Accuracy": [basevar_metrics["Accuracy"], glimpse_metrics["Accuracy"]],
        "Recall": [basevar_metrics["Recall"], glimpse_metrics["Recall"]],
        "F1 Score": [basevar_metrics["F1 Score"], glimpse_metrics["F1 Score"]],
    }

    df_stats = pd.DataFrame(stats_data)
    return df_stats


def generate_summary_statistics(df, output_dir):
    """
    Tạo thống kê cho toàn bộ DataFrame và biến thể hiếm (Rare Variants), trả về DataFrame.
    """
    try:
        logger.info("Generating summary statistics....")

        # Lọc biến thể hiếm
        rare_variants_df = df[df["INFO"].apply(is_rare_variant)]

        # Tính thống kê tổng quan và riêng cho biến thể hiếm
        total_stats_df = calculate_statistics(df, label="Total")
        rare_stats_df = calculate_statistics(rare_variants_df, label="Rare")

        # Gộp kết quả lại thành một DataFrame
        summary_df = pd.concat([total_stats_df, rare_stats_df], ignore_index=True)
        save_results_to_csv("sumary", summary_df, output_dir)

        logger.info("Summary statistics generated successfully")
        return summary_df
    except Exception as e:
        logger.error(f"Error generating summary statistics: {e}")
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

            df = compare_variants(ground_truth_path, basevar_path, glimpse_path, output_dir)
            generate_summary_statistics(df, output_dir)

        else:
            child, mom = sample_name.split("_")

            output_dir_mom = statistic_nipt_outdir(fq, chromosome, "mom")
            df = compare_variants(ground_truth_vcf(mom, chromosome), basevar_vcf(fq, chromosome), glimpse_vcf(fq, chromosome), output_dir_mom)
            generate_summary_statistics(df, output_dir_mom)

            output_dir_child = statistic_nipt_outdir(fq, chromosome, "child")
            df = compare_variants(ground_truth_vcf(child, chromosome), basevar_vcf(fq, chromosome), glimpse_vcf(fq, chromosome), output_dir_child)
            generate_summary_statistics(df, output_dir_child)

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
