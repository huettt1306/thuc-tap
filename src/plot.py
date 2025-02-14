import os
import pandas as pd
import matplotlib.pyplot as plt
from helper.config import PATHS, PARAMETERS, TRIO_DATA
from helper.file_utils import save_results_to_csv
from helper.path_define import statistic_outdir, fastq_single_path, fastq_nipt_path


def read_and_process_single_samples(chromosome, outdir):
    """
    Đọc và xử lý tất cả các mẫu.
    """
    print(f"Processing single sample....")
    for coverage in PARAMETERS["coverage"]:
        df_list = []
        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            for trio_name, trio_info in TRIO_DATA.items():
                for role, name in trio_info.items():

                    fq_path = os.path.join(fastq_single_path(name, coverage, index), f"{name}.fastq.gz")
                    sample_path = os.path.join(statistic_outdir(fq_path, chromosome), "summary.csv")
                    
                    if not os.path.exists(sample_path):
                        continue

                    df = pd.read_csv(sample_path)

                    # Đảm bảo các cột số được xử lý đúng kiểu
                    numeric_cols = df.columns.difference(['AF (%)'])  # Loại trừ cột 'AF (%)'
                    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

                    df["AF (%)"] = df["AF (%)"].str.replace("%", "").astype(float)
                    df_rare = df[(df["AF (%)"] > 0) & (df["AF (%)"] < 50)]

                    for col in df.columns:
                        df[f'Total {col}'] = df[col][::-1].cumsum()[::-1]
                        df[f'Rare {col}'] = df_rare[col][::-1].cumsum()[::-1]

                    category_list = ["AF (%)"]
                    for type in ["Total GT", "Total ALT", "Rare GT"]:
                        df[f"{type} correct ratio"] = df[f'{type} Glimpse True'] / df[f"{type} Glimpse"]
                        df[f"{type} found ratio"] = df[f'{type} Glimpse True'] / df[f"{type} Truth"]
                        df[f"{type} missing ratio"] = df[f'{type} Truth not found'] / df[f"{type} Truth"]
                        category_list.append(f"{type} correct ratio")
                        category_list.append(f"{type} found ratio")
                        category_list.append(f"{type} missing ratio")
                        
                    df_list.append(df)

        df_mean = pd.concat([df[category_list] for df in df_list]).groupby(level=0).mean()
        save_results_to_csv(os.path.join(outdir, f"summary_{coverage}x.csv"), df_mean)

        for type in ["Total GT", "Total ALT", "Rare GT"]:
            plot_mean_data(df_mean, os.path.join(outdir, f"{type}_{coverage}x.png"), type)

def read_and_process_nipt_samples(chromosome, outdir):
    """
    Đọc và xử lý tất cả các mẫu NIPT từ các thư mục trong root_dir.
    """
    print("Ploting nipt data...")
    for ff in PARAMETERS["ff"]:
        for coverage in PARAMETERS["coverage"]:
            df_list = []
            for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
                for trio_name, trio_info in TRIO_DATA.items():

                    child = trio_info["child"]
                    mother = trio_info["mother"]
                    father = trio_info["father"]

                    fq_path = os.path.join(fastq_nipt_path(child, mother, father, coverage, ff, index), f"{child}_{mother}_{father}.fastq.gz")
                    sample_path = os.path.join(statistic_outdir(fq_path, chromosome), "summary.csv")
                    
                    if not os.path.exists(sample_path):
                        continue

                    df = pd.read_csv(sample_path)

                    # Đảm bảo các cột số được xử lý đúng kiểu
                    numeric_cols = df.columns.difference(['AF (%)']) 
                    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

                    df["AF (%)"] = df["AF (%)"].str.replace("%", "").astype(float)
                    df_rare = df[(df["AF (%)"] > 0) & (df["AF (%)"] < 50)]

                    for col in df.columns:
                        df[f'Total {col}'] = df[col][::-1].cumsum()[::-1]
                        df[f'Rare {col}'] = df_rare[col][::-1].cumsum()[::-1]

                    category_list = ["AF (%)"]
                    for type in ["Total GT", "Total ALT", "Rare GT"]:
                        for truth in ["Mother", "Child"]:
                            df[f"{type} {truth} correct ratio"] = df[f'{type} Glimpse same as {truth}'] / df[f"{type} Glimpse"]
                            df[f"{type} {truth} found ratio"] = df[f'{type} Glimpse same as {truth}'] / df[f"{type} {truth}"]
                            df[f"{type} {truth} missing ratio"] = df[f'{type} {truth} not found'] / df[f"{type} {truth}"]
                            category_list.append(f"{type} {truth} correct ratio")
                            category_list.append(f"{type} {truth} found ratio")
                            category_list.append(f"{type} {truth} missing ratio")

                    df_list.append(df)

            df_mean = pd.concat([df[category_list] for df in df_list]).groupby(level=0).mean()
            save_results_to_csv(os.path.join(outdir, f"summary_{coverage}x_ff{ff}.csv"), df_mean)

            for type in ["Total GT", "Total ALT", "Rare GT"]:
                plot_mean_data(df_mean, os.path.join(outdir, f"{type}_{coverage}_x_ff{ff}.png"), type)


def plot_mean_data(df, outdir, category):
    x = df["AF (%)"]
    plt.figure(figsize=(20, 12))

    for column in df.columns:
        if column.startswith(category):
            plt.plot(x, df[column], marker='o', linestyle='-', label=column)

    plt.title(f"Biểu đồ Tỉ lệ {category} theo AF (%)", fontsize=14)
    plt.xlabel("AF (%)", fontsize=12)
    plt.ylabel("Giá trị", fontsize=12)
    plt.legend(title="Các cột", fontsize=10)
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(outdir, dpi=300)
    plt.close()

    print(f"Biểu đồ đã được lưu tại: {outdir}")


if __name__ == "__main__":
    for chromosome in PARAMETERS["chrs"]:
        
        output_dir = os.path.join(PATHS["plot_directory"], chromosome)
        os.makedirs(output_dir, exist_ok=True)

        read_and_process_single_samples(chromosome, "summary", output_dir)
        read_and_process_nipt_samples(chromosome, "summary", output_dir)
        print(f"Plots saved to {output_dir}")










































































