import pandas as pd
import os
from helper.GT import get_af_gt, get_af_gt_true, get_af_gt_false, get_af_gt_not_given
from helper.ALT import get_af_alt, get_af_alt_true, get_af_alt_false, get_af_alt_not_given
from helper.file_utils import save_results_to_csv, process_vcf
from helper.config import PATHS, PARAMETERS
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "single_statistic_pipeline.log"))


def compare_single_variants(ground_truth_path, basevar_path, glimpse_path, output_file):
    """
    So sánh biến thể giữa các phương pháp và lưu kết quả ra file CSV.
    """
    logger.info(f"Comparing variants for paths: {ground_truth_path}, {basevar_path}, {glimpse_path}")

    try:
        # Kiểm tra nếu file đã tồn tại, đọc DataFrame từ file
        if os.path.exists(output_file):
            logger.info(f"Output file {output_file} already exists. Loading existing results.")
            merged_df = pd.read_csv(output_file)
            return merged_df
        
        ground_truth_df = process_vcf(ground_truth_path, "Truth")
        basevar_df = process_vcf(basevar_path, "BaseVar")
        glimpse_df = process_vcf(glimpse_path, "Glimpse")

        merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")
        merged_df = pd.merge(ground_truth_df, merged_df, on=["CHROM", "POS", "REF", "ALT"], how="outer")

        merged_df.drop(columns=['AF_BaseVar'], inplace=True)
        merged_df.drop(columns=['AF_Glimpse'], inplace=True)
        merged_df.rename(columns={'AF_Truth': 'AF'}, inplace=True)

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
    if af < 0:
        return stats
    
    if af not in stats:
        stats[af] = {
            "GT Truth": 0,
            "GT Glimpse": 0,
            "GT Glimpse True": 0,
            "GT Glimpse False": 0,
            "GT Truth not found": 0,

            "ALT Truth": 0,
            "ALT Basevar": 0,
            "ALT Glimpse": 0,
            "ALT Glimpse True": 0,
            "ALT Glimpse False": 0,
            "ALT Truth not found": 0,
        }

    # Cập nhật giá trị của trường 'field'
    stats[af][field] += 1
    return stats


def calculate_af_single_statistics(df):
    """
    Tính toán thống kê cho từng giá trị AF và trả về dưới dạng dictionary.
    Mỗi entry trong dictionary sẽ chứa các thống kê cho một giá trị AF.
    """
    stats = {}

    for _, row in df.iterrows():
        stats = update_stats(stats, get_af_gt(row, "Truth"), "GT Truth")
        stats = update_stats(stats, get_af_gt(row, "Glimpse"), "GT Glimpse")
        stats = update_stats(stats, get_af_gt_true(row, "Glimpse", "Truth"), "GT Glimpse True")
        stats = update_stats(stats, get_af_gt_false(row, "Glimpse", "Truth"), "GT Glimpse False")
        stats = update_stats(stats, get_af_gt_not_given(row, "Truth", "Glimpse"), "GT Truth not found")

        stats = update_stats(stats, get_af_alt(row, "Truth"), "ALT Truth")
        stats = update_stats(stats, get_af_alt(row, "BaseVar"), "ALT Basevar")
        stats = update_stats(stats, get_af_alt(row, "Glimpse"), "ALT Glimpse")
        stats = update_stats(stats, get_af_alt_true(row, "Glimpse", "Truth"), "ALT Glimpse True")
        stats = update_stats(stats, get_af_alt_false(row, "Glimpse", "Truth"), "ALT Glimpse False")
        stats = update_stats(stats, get_af_alt_not_given(row, "Truth", "Glimpse"), "ALT Truth not found")

    return stats







