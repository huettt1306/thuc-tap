import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score, recall_score, f1_score
from helper.config import PATHS, PARAMETERS, TRIO_DATA
from helper.path_define import statistic_outdir, fastq_single_path, fastq_nipt_path

# Các hàm hỗ trợ
def is_rare(info_str):
    """
    Kiểm tra nếu biến thể là hiếm dựa trên giá trị AF_ESA trong chuỗi INFO.
    """
    try:
        info_parts = dict(item.split('=') for item in info_str.split(';') if '=' in item)
        af_esa = float(info_parts.get('AF_ESA', 1))  # Mặc định là 1 nếu không có AF_ESA
        return af_esa < 0.001
    except Exception:
        return False

def calculate_metrics(y_true, y_pred, prefix):
    """
    Tính toán Precision, Recall, F1 và trả về kết quả dưới dạng dict.
    """
    return {
        f'precision_{prefix}': precision_score(y_true, y_pred),
        f'recall_{prefix}': recall_score(y_true, y_pred),
        f'f1_{prefix}': f1_score(y_true, y_pred),
    }

def compute_general_statistics(df):
    """
    Tính toán các thông số Precision, Recall, F1 cho tất cả các biến thể.
    """
    stats = {}
    stats.update(calculate_metrics(df['truth'], df['glimpse'], 'glimpse'))
    stats.update(calculate_metrics(df['truth'], df['basevar'], 'basevar'))

    return stats

def compute_rare_statistics(df):
    """
    Tính toán các thông số Precision, Recall, F1 cho các biến thể hiếm.
    """
    # Lọc biến thể hiếm
    rare_df = df[df['INFO'].apply(is_rare)]
    if rare_df.empty:
        return {
            'precision_glimpse_rare': None,
            'recall_glimpse_rare': None,
            'f1_glimpse_rare': None,
            'precision_basevar_rare': None,
            'recall_basevar_rare': None,
            'f1_basevar_rare': None,
        }

    stats = {}
    stats.update(calculate_metrics(rare_df['truth'], rare_df['glimpse'], 'glimpse_rare'))
    stats.update(calculate_metrics(rare_df['truth'], rare_df['basevar'], 'basevar_rare'))

    return stats

def compute_statistics(df):
    """
    Tính toán các thông số Precision, Recall, F1 cho cả tổng thể và biến thể hiếm.
    """

    general_stats = compute_general_statistics(df)
    rare_stats = compute_rare_statistics(df)

    # Kết hợp kết quả
    return {**general_stats, **rare_stats}

# Hàm đọc và xử lý dữ liệu
def read_and_process_all_samples(chromosome):
    """
    Đọc và xử lý tất cả các mẫu từ các thư mục trong root_dir.
    """
    stats_list = []
    for coverage in PARAMETERS["coverage"]:
        combined_df = pd.DataFrame()
        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            for trio_name, trio_info in TRIO_DATA.items():
                for role, name in trio_info.items():
                    sample_path = os.path.join(statistic_outdir(fastq_single_path(name, coverage, index), chromosome), "variants.csv")
                    if not sample_path:
                        continue
                    sample_df = pd.read_csv(sample_path)
                    combined_df = pd.concat([combined_df, sample_df], ignore_index=True)

        stats = compute_statistics(combined_df)
        stats['cov'] = coverage
        stats_list.append(stats) 

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

if __name__ == "__main__":
    for chromosome in PARAMETERS["chrs"]:
        
        # Đường dẫn để lưu các biểu đồ
        output_dir = os.path.join(PATHS["plot_directory"], chromosome)
        os.makedirs(output_dir, exist_ok=True)

        # Đọc và xử lý dữ liệu từ tất cả các mẫu
        stats_df = read_and_process_all_samples(chromosome)

        # Vẽ biểu đồ
        plot_statistics(stats_df, output_dir)

        print(f"Plots saved to {output_dir}")
