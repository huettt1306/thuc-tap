import pandas as pd
import os
from cyvcf2 import VCF
from helper.converter import convert_genotype
from helper.file_utils import save_results_to_csv
from helper.path_define import ground_truth_vcf, statistic_variants, statistic_summary, glimpse_vcf, basevar_vcf, samid
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
                f"ALT_{method_name}": str(record.ALT[0]) if record.ALT else None,
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

        merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF"], how="outer")
        merged_df = pd.merge(ground_truth_df, merged_df, on=["CHROM", "POS", "REF"], how="outer")

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

        merged_df = pd.merge(glimpse_df, basevar_df, on=["CHROM", "POS", "REF"], how="outer")
        merged_df = pd.merge(child_df, merged_df, on=["CHROM", "POS", "REF"], how="outer")
        merged_df = pd.merge(mother_df, merged_df, on=["CHROM", "POS", "REF"], how="outer")
        merged_df = pd.merge(father_df, merged_df, on=["CHROM", "POS", "REF"], how="outer")
        
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


def count_row_variants(row, method, truth):
    """
    Đếm số lượng biến thể theo dòng.
    """
    variants = 0
    variants_gt = 0

    if row[f"ALT_{method}"] == row[f"ALT_{truth}"]:
        variants += 1
        if row[f"GT_{method}"] == row[f"GT_{truth}"]:
            variants_gt += 1

    return variants, variants_gt



def calculate_af_statistics(df, af_percent):
    """
    Tính toán thống kê cho một giá trị AF cho trước và trả về dưới dạng dictionary.
    """
    stats = {
        "Total Variants": 0,
        "Truth Variants": 0,
        "BaseVar Variants": 0,
        "BaseVar Variants True": 0,
        "BaseVar Variants GT True": 0,
        "Glimpse Variants": 0,
        "Glimpse Variants True": 0,
        "Glimpse Variants GT True": 0,
    }

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        af_value = row["AF"] 
        # Kiểm tra xem AF có nằm trong khoảng từ 0% đến af_percent%
        if 0 < af_value <= af_percent / 100:
            stats["Total Variants"] += 1

            if row["BaseVar"]:
                stats["BaseVar Variants"] += 1

            if row["Glimpse"]:
                stats["Glimpse Variants"] += 1

            if row["Truth"]:
                stats["Truth Variants"] += 1

                # Tính toán cho BaseVar
                var_true, var_gt = count_row_variants(row, "BaseVar", "Truth")
                stats["BaseVar Variants True"] += var_true
                stats["BaseVar Variants GT True"] += var_gt

                # Tính toán cho Glimpse
                var_true, var_gt = count_row_variants(row, "Glimpse", "Truth")
                stats["Glimpse Variants True"] += var_true
                stats["Glimpse Variants GT True"] += var_gt

    return stats


def calculate_af_nipt_statistics(df, af_percent):
    """
    Tính toán thống kê cho một giá trị AF cho trước.
    """
    stats = {
        "Total Variants": 0,
        "Child Variants": 0,
        "Mother Variants": 0,
        "Father Variants": 0,
        "BaseVar Variants": 0,
        "BaseVar Variants Child": 0,
        "BaseVar Variants Mother": 0,
        "BaseVar Variants Father": 0,
        "Glimpse Variants": 0,
        "Glimpse Variants Child": 0,
        "Glimpse Variants Child GT": 0,
        "Glimpse Variants Mother": 0,
        "Glimpse Variants Mother GT": 0,
        "Glimpse Variants Father": 0,
        "Glimpse Variants Father GT": 0,
        "Priv Father Variants": 0,
        "Priv Glimpse Variants Father": 0,
        "Priv Mother Variants": 0,
        "Priv BaseVar Variants Mother": 0,
        "Priv Glimpse Variants Mother": 0,
        "Diff Glimpse Variants Parent GT": 0,
        "Diff Variants Parent GT": 0,
    }

    # Duyệt qua từng dòng trong DataFrame df
    for _, row in df.iterrows():
        af_value = row["AF"] 
        if 0 < af_value <= af_percent / 100:
            stats["Total Variants"] += 1
            if row["BaseVar"]:
                stats["BaseVar Variants"] += 1

            if row["Glimpse"]:
                stats["Glimpse Variants"] += 1

            if row["Child"]:
                stats["Child Variants"] += 1

                basevar_true, basevar_gt = count_row_variants(row, "BaseVar", "Child")
                stats["BaseVar Variants Child"] += basevar_true

                glimpse_true, glimpse_gt = count_row_variants(row, "Glimpse", "Child")
                stats["Glimpse Variants Child"] += glimpse_true
                stats["Glimpse Variants Child GT"] += glimpse_gt

            if row["Mother"]:
                stats["Mother Variants"] += 1

                basevar_true, basevar_gt = count_row_variants(row, "BaseVar", "Mother")
                stats["BaseVar Variants Mother"] += basevar_true

                glimpse_true, glimpse_gt = count_row_variants(row, "Glimpse", "Mother")
                stats["Glimpse Variants Mother"] += glimpse_true
                stats["Glimpse Variants Mother GT"] += glimpse_gt

                if not row["Father"]:
                    stats["Priv Mother Variants"] += 1
                    stats["Priv BaseVar Variants Mother"] += basevar_true
                    stats["Priv Glimpse Variants Mother"] += glimpse_true

                elif row["GT_Father"] != row["GT_Mother"]:
                    stats["Diff Variants Parent GT"] += 1
                    stats["Diff Glimpse Variants Parent GT"] += glimpse_gt

            if row["Father"]:
                stats["Father Variants"] += 1

                basevar_true, basevar_gt = count_row_variants(row, "BaseVar", "Father")
                stats["BaseVar Variants Father"] += basevar_true

                glimpse_true, glimpse_gt = count_row_variants(row, "Glimpse", "Father")
                stats["Glimpse Variants Father"] += glimpse_true
                stats["Glimpse Variants Father GT"] += glimpse_gt

                if not row["Mother"]:
                    stats["Priv Father Variants"] += 1
                    stats["Priv Glimpse Variants Father"] += glimpse_true

    return stats


def generate_summary_statistics(df, output_file, type="single"):
    """
    Tạo thống kê cho toàn bộ DataFrame theo từng AF.
    """
    try:
        logger.info("Generating summary statistics....")

        # Khởi tạo stats_data
        stats_data = {"AF (%)": []}

        # Danh sách các giá trị AF cần tính
        af_values = PARAMETERS["af"]

        # Duyệt qua các giá trị AF đã định nghĩa
        for af_percent in af_values:
            logger.info(f"Generating summary statistics for {af_percent}%....")
            if type == "nipt" :
                stats = calculate_af_nipt_statistics(df, af_percent)
            else :
                stats = calculate_af_statistics(df, af_percent)

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
            return
            df = compare_variants(ground_truth_vcf(sample_name, chromosome), basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            generate_summary_statistics(df, statistic_summary(fq, chromosome))

        else:
            child, mom, dad = sample_name.split("_")
            child_path = ground_truth_vcf(child, chromosome)
            mother_path = ground_truth_vcf(mom, chromosome)
            father_path = ground_truth_vcf(dad, chromosome)

            df = compare_nipt_variants(child_path, mother_path, father_path, basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            generate_summary_statistics(df, statistic_summary(fq, chromosome), "nipt")
            
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
