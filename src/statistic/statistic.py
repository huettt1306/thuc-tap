import pandas as pd
import os
from cyvcf2 import VCF
from helper.variant import convert_genotype, count_variant, count_same_gt, count_priv_gt, count_priv_true_gt, count_same_true_gt, count_same_false_gt
from helper.file_utils import save_results_to_csv
from helper.path_define import ground_truth_vcf, statistic_variants, statistic_summary, glimpse_vcf, basevar_vcf, samid, statistic_rare_summary
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
            af_value = dict(record.INFO).get('AF', None)
            variants.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT": str(record.ALT[0]) if record.ALT else None,
                f"AF_{method_name}": af_value[0] if isinstance(af_value, tuple) else af_value if isinstance(af_value, float) else None,
                f"GT_{method_name}": convert_genotype(record.genotypes[0]),
                method_name: True
            })
    except Exception as e:
        logger.error(f"Error processing VCF file {vcf_path}: {e}")
        raise

    logger.info(f"Finished processing VCF file: {vcf_path}")
    return pd.DataFrame(variants)

def compare_variants(ground_truth_path, basevar_path, glimpse_path, output_file):
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


def calculate_af_statistics(df, af_percent, af_rare):
    """
    Tính toán thống kê cho một giá trị AF cho trước và trả về dưới dạng dictionary.
    """
    stats = {
        "Truth Variants": 0,
        "BaseVar Variants": 0,
        "BaseVar True Variants": 0,
        "Glimpse Variants": 0,
        "Glimpse True Variants": 0,
    }

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        af = af_percent / 100
        stats["Truth Variants"] += count_variant(row, "Truth", af, af_rare)
        stats["BaseVar Variants"] += count_variant(row, "BaseVar", af, af_rare)
        stats["Glimpse Variants"] += count_variant(row, "Glimpse", af, af_rare)
        stats["BaseVar True Variants"] += count_same_gt(row, "BaseVar", "Truth", af, af_rare)
        stats["Glimpse True Variants"] += count_same_gt(row, "Glimpse", "Truth", af, af_rare)

    return stats


def calculate_af_nipt_statistics(df, af_percent, af_rare):
    """
    Tính toán thống kê cho một giá trị AF cho trước.
    """
    stats = {
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

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        af = af_percent / 100
        stats["Child Variants"] += count_variant(row, "Child", af, af_rare)
        stats["Mother Variants"] += count_variant(row, "Mother", af, af_rare)
        stats["Father Variants"] += count_variant(row, "Father", af, af_rare)
        stats["BaseVar Variants"] += count_variant(row, "BaseVar", af, af_rare)
        stats["Glimpse Variants"] += count_variant(row, "Glimpse", af, af_rare)

        stats["Child Variants different from Mother"] += count_priv_gt(row, "Child", "Mother", af, af_rare)
        stats["Father Variants different from Mother"] += count_priv_gt(row, "Father", "Mother", af, af_rare)
        stats["Child Variants same as Mother"] += count_same_gt(row, "Child", "Mother", af, af_rare)
        stats["Father Variants same as Mother"] += count_same_gt(row, "Father", "Mother", af, af_rare)

        stats["Glimpse Variants same as Child"] += count_same_gt(row, "Glimpse", "Child", af, af_rare)
        stats["Glimpse Variants same as Mother"] += count_same_gt(row, "Glimpse", "Mother", af, af_rare)

        stats["Glimpse Variants same as Child but different from Mother"] += count_priv_true_gt(row, "Glimpse", "Child", "Mother", af, af_rare)
        stats["Glimpse Variants same as Mother but different from Child"] += count_priv_true_gt(row, "Glimpse", "Mother", "Child", af, af_rare)
        stats["Glimpse Variants same as Child and Mother"] += count_same_true_gt(row, "Glimpse", "Child", "Mother", af, af_rare)
        stats["Glimpse Variants Different from Child and Mother"] += count_same_false_gt(row, "Glimpse", "Child", "Mother", af, af_rare)
        
    return stats


def generate_summary_statistics(df, output_file, af_rare, type="single"):
    """
    Tạo thống kê cho toàn bộ DataFrame theo từng AF.
    """
    try:
        logger.info("Generating summary statistics....")

        # Khởi tạo stats_data
        stats_data = {"AF (%)": []}

        # Duyệt qua các giá trị AF đã định nghĩa
        for af_percent in PARAMETERS["af"]:
            logger.info(f"Generating summary statistics for {af_percent}%....")
            if type == "nipt" :
                stats = calculate_af_nipt_statistics(df, af_percent, af_rare)
            else :
                stats = calculate_af_statistics(df, af_percent, af_rare)

            # Cập nhật giá trị AF
            stats_data["AF (%)"].append(f"{af_percent}%")

            # Cập nhật các giá trị thống kê vào stats_data
            for key, value in stats.items():
                if key not in stats_data:
                    stats_data[key] = []
                stats_data[key].append(value)

        # Tạo DataFrame từ stats_data
        summary_df = pd.DataFrame(stats_data)

        # Lưu kết quả vào file CSV
        save_results_to_csv(output_file, summary_df)

        logger.info("Summary statistics generated successfully")
        return summary_df

    except Exception as e:
        logger.error(f"Error generating summary statistics: {e}")
        raise


def statistic(fq, chromosome):
    sample_name = samid(fq)

    try:
        logger.info(f"Starting statistics for sample {sample_name} on chromosome {chromosome}")
        
        basevar_path = basevar_vcf(fq, chromosome)
        glimpse_path = glimpse_vcf(fq, chromosome)

        if "_" not in sample_name:
            df = compare_variants(ground_truth_vcf(sample_name, chromosome), basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            generate_summary_statistics(df, statistic_summary(fq, chromosome), 1)
            generate_summary_statistics(df, statistic_rare_summary(fq, chromosome), PARAMETERS["af_rare"])
            
        else:
            child, mom, dad = sample_name.split("_")
            child_path = ground_truth_vcf(child, chromosome)
            mother_path = ground_truth_vcf(mom, chromosome)
            father_path = ground_truth_vcf(dad, chromosome)

            df = compare_nipt_variants(child_path, mother_path, father_path, basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            generate_summary_statistics(df, statistic_summary(fq, chromosome), 1, "nipt")
            generate_summary_statistics(df, statistic_rare_summary(fq, chromosome), PARAMETERS["af_rare"], "nipt")
            
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
