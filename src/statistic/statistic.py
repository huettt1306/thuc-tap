import pandas as pd
import os
from helper.file_utils import save_results_to_csv
from helper.path_define import ground_truth_vcf, statistic_variants, statistic_summary, glimpse_vcf, basevar_vcf, samid
from helper.config import PATHS, PARAMETERS
from helper.logger import setup_logger
from statistic.single_stats import compare_single_variants, calculate_af_single_statistics
from statistic.nipt_stats import compare_nipt_variants, calculate_af_nipt_statistics

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "statistic_pipeline.log"))


def process_dataframe(df):
    # Chuyển tất cả các cột sang kiểu số, nếu lỗi -> thay bằng 0
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Sắp xếp theo AF (%) tăng dần
    df = df.sort_values(by="AF (%)", ascending=True)

    # Tính tổng tích lũy cho tất cả các cột
    for col in df.columns:
        if col != "AF (%)":
            df[f'Total {col}'] = df[col][::-1].cumsum()[::-1]

    return df


def generate_summary_statistics(df, output_file, type="single"):
    """
    Tạo thống kê cho toàn bộ DataFrame theo từng AF.
    """
    try:
        logger.info("Generating summary statistics....")

        # Khởi tạo stats_data
        stats_data = {"AF (%)": []}

        # Tính toán thống kê cho từng loại AF
        if type == "nipt":
            stats = calculate_af_nipt_statistics(df)
        else:
            stats = calculate_af_single_statistics(df)

        # Duyệt qua các AF trong kết quả thống kê đã tính toán
        for af_percent, stat_values in stats.items():
            # Cập nhật giá trị AF vào stats_data
            stats_data["AF (%)"].append(int(af_percent))

            # Cập nhật các giá trị thống kê vào stats_data
            for key, value in stat_values.items():
                if key not in stats_data:
                    stats_data[key] = []
                stats_data[key].append(value)

        # Tạo DataFrame từ stats_data
        summary_df = process_dataframe(pd.DataFrame(stats_data))

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
            df = compare_single_variants(ground_truth_vcf(sample_name, chromosome), basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            return generate_summary_statistics(df, statistic_summary(fq, chromosome))
            
        else:
            child, mom, dad = sample_name.split("_")
            child_path = ground_truth_vcf(child, chromosome)
            mother_path = ground_truth_vcf(mom, chromosome)
            father_path = ground_truth_vcf(dad, chromosome)

            df = compare_nipt_variants(child_path, mother_path, father_path, basevar_path, glimpse_path, statistic_variants(fq, chromosome))
            return generate_summary_statistics(df, statistic_summary(fq, chromosome), "nipt")
            
    except Exception as e:
        logger.error(f"Error in statistic function for sample {sample_name} and chromosome {chromosome}: {e}")
        raise


def run_statistic(fq):
    try:
        logger.info(f"Starting full statistics pipeline for sample {fq}")

        all_data = []  # Dùng list để lưu các DataFrame tạm thời

        for chromosome in PARAMETERS["chrs"]:
            df = statistic(fq, chromosome)

            if df is not None and not df.empty:
                all_data.append(df)  # Thêm vào danh sách, chưa merge ngay

        if all_data:  # Nếu có dữ liệu, thực hiện gộp và cộng tổng
            all_df = pd.concat(all_data).groupby("AF", as_index=False).sum()
        else:
            all_df = pd.DataFrame()  # Nếu không có dữ liệu, trả về DataFrame rỗng

        save_results_to_csv(statistic_summary(fq), all_df)
        logger.info(f"Completed full statistics pipeline for sample {fq}")

        return all_df
    except Exception as e:
        logger.error(f"Error in run_statistic pipeline for sample {fq}: {e}")
        raise