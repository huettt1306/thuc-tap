import os
import pandas as pd
import matplotlib.pyplot as plt
from helper.config import PATHS, PARAMETERS, TRIO_DATA
from helper.path_define import statistic_outdir, statistic_nipt_outdir, fastq_single_path, fastq_nipt_path

# Hàm đọc và xử lý dữ liệu
def read_and_process_single_samples(chromosome):
    """
    Đọc và xử lý tất cả các mẫu từ các thư mục trong root_dir.
    """
    stats_list = []
    for coverage in PARAMETERS["coverage"]:
        combined_df = pd.DataFrame()
        sample_dfs = []

        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            for trio_name, trio_info in TRIO_DATA.items():
                for role, name in trio_info.items():
                    sample_path = os.path.join(statistic_outdir(fastq_single_path(name, coverage, index), chromosome), "summary.csv")
                    if not os.path.exists(sample_path):
                        continue

# làm gì đó

    return pd.DataFrame(stats_list)

def read_and_process_nipt_samples(chromosome):
    """
    Đọc và xử lý tất cả các mẫu NIPT từ các thư mục trong root_dir.
    """
    stats_list = []
    for coverage in PARAMETERS["coverage"]:
        for ff in PARAMETERS["ff"]:
            child_dfs = []  # Danh sách các DataFrame cho con
            mother_dfs = []  # Danh sách các DataFrame cho mẹ

            for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
                for trio_name, trio_info in TRIO_DATA.items():
                    # Đường dẫn đến dữ liệu của con và mẹ
                    sample = fastq_nipt_path(trio_info["child"], trio_info["mother"], coverage, ff, index)
                    child_path = os.path.join(statistic_nipt_outdir(sample, chromosome, "child"), "variants.csv")
                    mother_path = os.path.join(statistic_nipt_outdir(sample, chromosome, "mother"), "variants.csv")

                    # Đọc dữ liệu con
                    if os.path.exists(child_path):
                        child_df = pd.read_csv(child_path)
                        child_dfs.append(child_df)

                    # Đọc dữ liệu mẹ
                    if os.path.exists(mother_path):
                        mother_df = pd.read_csv(mother_path)
                        mother_dfs.append(mother_df)

# làm gì đó
            
    return pd.DataFrame(stats_list)


# Hàm vẽ biểu đồ
def plot_statistics(stats_df, output_dir):
    """
    Vẽ các biểu đồ Precision, Recall, F1 cho Glimpse và Basevar.
    """
    metrics = ['precision', 'recall', 'f1']
    prefixes = ['', '_rare']
    methods = ['glimpse', 'basevar']

    for metric in metrics:
        for prefix in prefixes:
            plt.figure(figsize=(10, 6))
            for method in methods:
                column_name = f'{metric}_{method}{prefix}'
                if column_name in stats_df.columns:
                    plt.plot(stats_df['cov'], stats_df[column_name], marker='o', label=method.capitalize())
            
            # Thiết lập biểu đồ
            plt.title(f'{metric.capitalize()} Comparison{" (Rare Variants)" if prefix else ""}')
            plt.xlabel('Coverage')
            plt.ylabel(metric.capitalize())
            plt.legend()
            plt.grid(True)
            plt.tight_layout()

            # Lưu biểu đồ
            plot_name = f'{metric}{"_rare" if prefix else ""}.png'
            output_path = os.path.join(output_dir, plot_name)
            plt.savefig(output_path)
            plt.close()

## viết lại hàm plot_statistics_nipt tương ứng. Tôi cần 3 biểu đồ là Bubble Chart, Line Plot và Heatmap thể hiện tương quan giữa cov, ff với 3 chỉ số thống kê. Vẽ các biểu đồ riêng khi so sánh với ground-truth mẹ và con, tương ứng với code read_and_process_nipt_samples bên trên của bạn

if __name__ == "__main__":
    for chromosome in PARAMETERS["chrs"]:
        
        # Đường dẫn để lưu các biểu đồ
        output_dir = os.path.join(PATHS["plot_directory"], chromosome)
        os.makedirs(output_dir, exist_ok=True)

        # Đọc và xử lý dữ liệu từ tất cả các mẫu
        stats_df = read_and_process_single_samples(chromosome)

        # Vẽ biểu đồ
        plot_statistics(stats_df, output_dir)

        print(f"Plots saved to {output_dir}")

