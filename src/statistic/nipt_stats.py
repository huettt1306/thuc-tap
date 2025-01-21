import pandas as pd
import os
from cyvcf2 import VCF
from helper.variant import af_variant, af_same_gt, af_priv_gt, af_priv_true_gt, af_same_false_gt, af_same_true_gt
from helper.file_utils import save_results_to_csv, process_vcf
from helper.config import PATHS, PARAMETERS
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "nipt_statistic_pipeline.log"))


def compare_nipt_variants(child_path, mother_path, father_path, basevar_path, glimpse_path, output_file):
    """
    So sánh biến thể giữa các phương pháp và lưu kết quả ra file CSV.
    """
    logger.info(f"Comparing variants for paths: {basevar_path}, {glimpse_path}")

    try:
        if os.path.exists(output_file):
            logger.info(f"Output file {output_file} already exists. Loading existing results.")
            merged_df = pd.read_csv(output_file)
            return merged_df
        
        child_df = process_vcf(child_path, "Child")
        mother_df = process_vcf(mother_path, "Mother")
        father_df = process_vcf(father_path, "Father")
        basevar_df = process_vcf(basevar_path, "BaseVar")
        glimpse_df = process_vcf(glimpse_path, "Glimpse")

        merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df = pd.merge(child_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df = pd.merge(mother_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df = pd.merge(father_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        
        merged_df['AF'] = merged_df['AF_Mother'].combine_first(merged_df['AF_Father'])
        merged_df.drop(columns=['AF_BaseVar'], inplace=True)
        merged_df.drop(columns=['AF_Glimpse'], inplace=True)
        merged_df.drop(columns=['AF_Mother'], inplace=True)
        merged_df.drop(columns=['AF_Father'], inplace=True)
        merged_df.drop(columns=['AF_Child'], inplace=True)
        
        merged_df["AF"].fillna(-1.0, inplace=True)
        merged_df.fillna(False, inplace=True)

        save_results_to_csv(output_file, merged_df)
        logger.info(f"Detailed variant comparison and summary statistics saved to {output_file}")

        return merged_df
    except Exception as e:
        logger.error(f"Error comparing variants: {e}")
        raise


def update_stats(stats, af, field):
    """
    Kiểm tra xem af đã có trong stats chưa, nếu chưa thì thêm mới, nếu có thì cập nhật trường 'field'.
    """
    if af not in stats:
        # Nếu chưa có, khởi tạo thống kê cho af
        stats[af] = {
            "Child Variants": 0,
            "Mother Variants": 0,
            "Father Variants": 0,
            "BaseVar Variants": 0,
            "Glimpse Variants": 0,

            "Child Variants different from Mother": 0,
            "Father Variants different from Mother": 0,
            "Child Variants same as Mother": 0,
            "Father Variants same as Mother": 0,

            "Glimpse Variants same as Child": 0,
            "Glimpse Variants same as Mother": 0,
            "Glimpse Variants same as Child but different from Mother": 0,
            "Glimpse Variants same as Mother but different from Child": 0,
            "Glimpse Variants same as Child and Mother": 0,
            "Glimpse Variants Different from Child and Mother": 0,
        }

    # Cập nhật giá trị của trường 'field'
    stats[af][field] += 1
    return stats


def calculate_af_nipt_statistics(df):
    """
    Tính toán thống kê cho một giá trị AF cho trước.
    """
    stats = {}

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        stats = update_stats(stats, af_variant(row, "Child"), "Child Variants")
        stats = update_stats(stats, af_variant(row, "Mother"), "Mother Variants")
        stats = update_stats(stats, af_variant(row, "Father"), "Father Variants")
        stats = update_stats(stats, af_variant(row, "BaseVar"), "BaseVar Variants")
        stats = update_stats(stats, af_variant(row, "Glimpse"), "Glimpse Variants")

        stats = update_stats(stats, af_priv_gt(row, "Child", "Mother"), "Child Variants different from Mother")
        stats = update_stats(stats, af_priv_gt(row, "Father", "Mother"), "Father Variants different from Mother")
        stats = update_stats(stats, af_same_gt(row, "Child", "Mother"), "Child Variants same as Mother")
        stats = update_stats(stats, af_same_gt(row, "Father", "Mother"), "Father Variants same as Mother")

        stats = update_stats(stats, af_same_gt(row, "Glimpse", "Child"), "Glimpse Variants same as Child")
        stats = update_stats(stats, af_same_gt(row, "Glimpse", "Mother"), "Glimpse Variants same as Mother")

        stats = update_stats(stats, af_priv_true_gt(row, "Glimpse", "Child", "Mother"), "Glimpse Variants same as Child but different from Mother")
        stats = update_stats(stats, af_priv_true_gt(row, "Glimpse", "Mother", "Child"), "Glimpse Variants same as Mother but different from Child")
        stats = update_stats(stats, af_same_true_gt(row, "Glimpse", "Child", "Mother"), "Glimpse Variants same as Child and Mother")
        stats = update_stats(stats, af_same_false_gt(row, "Glimpse", "Child", "Mother"), "Glimpse Variants Different from Child and Mother")
        
    return stats







