import pandas as pd
import os
from cyvcf2 import VCF
from helper.variant import af_variant, af_same_gt
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
    if af not in stats:
        # Nếu chưa có, khởi tạo thống kê cho af
        stats[af] = {
            "Truth Variants": 0,
            "BaseVar Variants": 0,
            "BaseVar True Variants": 0,
            "Glimpse Variants": 0,
            "Glimpse True Variants": 0
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

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        # Cập nhật thống kê cho từng loại
        stats = update_stats(stats, af_variant(row, "Truth"), "Truth Variants")
        stats = update_stats(stats, af_variant(row, "BaseVar"), "BaseVar Variants")
        stats = update_stats(stats, af_variant(row, "Glimpse"), "Glimpse Variants")
        stats = update_stats(stats, af_same_gt(row, "BaseVar", "Truth"), "BaseVar True Variants")
        stats = update_stats(stats, af_same_gt(row, "Glimpse", "Truth"), "Glimpse True Variants")

    return stats







