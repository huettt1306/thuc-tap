import pandas as pd
import os
from helper.GT import get_af_gt, get_af_gt_true, get_af_gt_false, get_af_gt_not_given, get_af_gt_priv_true, get_af_gt_same_true, get_af_gt_same_false
from helper.ALT import get_af_alt, get_af_alt_true, get_af_alt_false, get_af_alt_not_given, get_af_alt_priv_true, get_af_alt_same_true, get_af_alt_same_false
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
    if af < 0:
        return stats
    
    if af not in stats:
        # Nếu chưa có, khởi tạo thống kê cho af
        stats[af] = {
            "GT Child": 0,
            "GT Mother": 0,
            "GT Father": 0,
            "GT Glimpse": 0,

            "GT Child diff from Mother": 0,
            "GT Father diff from Mother": 0,
            "GT Child same as Mother": 0,
            "GT Child same as Father": 0,
            "GT Father same as Mother": 0,

            "GT Glimpse same as Child": 0,
            "GT Glimpse same as Mother": 0,
            "GT Glimpse diff from Child": 0,
            "GT Glimpse diff from Mother": 0,

            "GT Glimpse same as Child but diff from Mother": 0,
            "GT Glimpse same as Mother but diff from Child": 0,
            "GT Glimpse same as Child and Mother": 0,
            "GT Glimpse diff from Child and Mother": 0,

            "GT Child not found": 0,
            "GT Mother not found": 0,
            "ALT Child not found": 0,
            "ALT Mother not found": 0,
            
            "ALT Child": 0,
            "ALT Mother": 0,
            "ALT Father": 0,
            "ALT Glimpse": 0,

            "ALT Child diff from Mother": 0,
            "ALT Father diff from Mother": 0,
            "ALT Child same as Mother": 0,
            "ALT Child same as Father": 0,
            "ALT Father same as Mother": 0,

            "ALT Glimpse same as Child": 0,
            "ALT Glimpse same as Mother": 0,
            "ALT Glimpse diff from Child": 0,
            "ALT Glimpse diff from Mother": 0,

            "ALT Glimpse same as Child but diff from Mother": 0,
            "ALT Glimpse same as Mother but diff from Child": 0,
            "ALT Glimpse same as Child and Mother": 0,
            "ALT Glimpse diff from Child and Mother": 0,
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
        stats = update_stats(stats, get_af_gt(row, "Child"), "GT Child")
        stats = update_stats(stats, get_af_gt(row, "Mother"), "GT Mother")
        stats = update_stats(stats, get_af_gt(row, "Father"), "GT Father")
        stats = update_stats(stats, get_af_gt(row, "Glimpse"), "GT Glimpse")

        stats = update_stats(stats, get_af_gt_false(row, "Child", "Mother"), "GT Child diff from Mother")
        stats = update_stats(stats, get_af_gt_false(row, "Father", "Mother"), "GT Father diff from Mother")
        stats = update_stats(stats, get_af_gt_true(row, "Child", "Mother"), "GT Child same as Mother")
        stats = update_stats(stats, get_af_gt_true(row, "Child", "Father"), "GT Child same as Father")
        stats = update_stats(stats, get_af_gt_true(row, "Father", "Mother"), "GT Father same as Mother")

        stats = update_stats(stats, get_af_gt_true(row, "Glimpse", "Child"), "GT Glimpse same as Child")
        stats = update_stats(stats, get_af_gt_true(row, "Glimpse", "Mother"), "GT Glimpse same as Mother")
        stats = update_stats(stats, get_af_gt_false(row, "Glimpse", "Child"), "GT Glimpse diff from Child")
        stats = update_stats(stats, get_af_gt_false(row, "Glimpse", "Mother"), "GT Glimpse diff from Mother")

        stats = update_stats(stats, get_af_gt_priv_true(row, "Glimpse", "Child", "Mother"), "GT Glimpse same as Child but diff from Mother")
        stats = update_stats(stats, get_af_gt_priv_true(row, "Glimpse", "Mother", "Child"), "GT Glimpse same as Mother but diff from Child")
        stats = update_stats(stats, get_af_gt_same_true(row, "Glimpse", "Child", "Mother"), "GT Glimpse same as Child and Mother")
        stats = update_stats(stats, get_af_gt_same_false(row, "Glimpse", "Child", "Mother"), "GT Glimpse diff from Child and Mother")
        
        stats = update_stats(stats, get_af_gt_not_given(row, "Child", "Glimpse"), "GT Child not found")
        stats = update_stats(stats, get_af_gt_not_given(row, "Mother", "Glimpse"), "GT Mother not found")
        stats = update_stats(stats, get_af_alt_not_given(row, "Child", "Glimpse"), "ALT Child not found")
        stats = update_stats(stats, get_af_alt_not_given(row, "Mother", "Glimpse"), "ALT Mother not found")
        
        stats = update_stats(stats, get_af_alt(row, "Child"), "ALT Child")
        stats = update_stats(stats, get_af_alt(row, "Mother"), "ALT Mother")
        stats = update_stats(stats, get_af_alt(row, "Father"), "ALT Father")
        stats = update_stats(stats, get_af_alt(row, "Glimpse"), "ALT Glimpse")

        stats = update_stats(stats, get_af_alt_false(row, "Child", "Mother"), "ALT Child diff from Mother")
        stats = update_stats(stats, get_af_alt_false(row, "Father", "Mother"), "ALT Father diff from Mother")
        stats = update_stats(stats, get_af_alt_true(row, "Child", "Mother"), "ALT Child same as Mother")
        stats = update_stats(stats, get_af_alt_true(row, "Child", "Father"), "ALT Child same as Father")
        stats = update_stats(stats, get_af_alt_true(row, "Father", "Mother"), "ALT Father same as Mother")

        stats = update_stats(stats, get_af_alt_true(row, "Glimpse", "Child"), "ALT Glimpse same as Child")
        stats = update_stats(stats, get_af_alt_true(row, "Glimpse", "Mother"), "ALT Glimpse same as Mother")
        stats = update_stats(stats, get_af_alt_false(row, "Glimpse", "Child"), "ALT Glimpse diff from Child")
        stats = update_stats(stats, get_af_alt_false(row, "Glimpse", "Mother"), "ALT Glimpse diff from Mother")

        stats = update_stats(stats, get_af_alt_priv_true(row, "Glimpse", "Child", "Mother"), "ALT Glimpse same as Child but diff from Mother")
        stats = update_stats(stats, get_af_alt_priv_true(row, "Glimpse", "Mother", "Child"), "ALT Glimpse same as Mother but diff from Child")
        stats = update_stats(stats, get_af_alt_same_true(row, "Glimpse", "Child", "Mother"), "ALT Glimpse same as Child and Mother")
        stats = update_stats(stats, get_af_alt_same_false(row, "Glimpse", "Child", "Mother"), "ALT Glimpse diff from Child and Mother")

    return stats







