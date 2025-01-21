import os
import pandas as pd
import matplotlib.pyplot as plt
from helper.config import PATHS, PARAMETERS, TRIO_DATA
from helper.file_utils import save_results_to_csv
from helper.path_define import statistic_outdir, fastq_single_path, fastq_nipt_path


# Hàm đọc và xử lý dữ liệu
def read_and_process_single_samples(chromosome, category, outdir):
    """
    Đọc và xử lý tất cả các mẫu từ các thư mục trong root_dir.
    """
    result = pd.DataFrame()
    for coverage in PARAMETERS["coverage"]:
        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            for trio_name, trio_info in TRIO_DATA.items():
                for role, name in trio_info.items():
                    # Đường dẫn tới tệp tin
                    fq_path = os.path.join(fastq_single_path(name, coverage, index), f"{name}.fastq.gz")
                    sample_path = os.path.join(statistic_outdir(fq_path, chromosome), f"{category}.csv")
                    
                    # Kiểm tra sự tồn tại của file
                    if not os.path.exists(sample_path):
                        continue

                    # Đọc file CSV
                    df = pd.read_csv(sample_path)

                    # Đảm bảo các cột số được xử lý đúng kiểu
                    numeric_cols = df.columns.difference(['AF (%)'])  # Loại trừ cột 'AF (%)'
                    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

                    # Nếu `result` trống, khởi tạo bằng `df`
                    if result.empty:
                        result = df
                    else:
                        # Cộng dồn dữ liệu
                        result = result.set_index('AF (%)').add(df.set_index('AF (%)'), fill_value=0).reset_index()

    result["Glimpse Accuracy"] = df["Glimpse True Variants"] / df["Glimpse Variants"]
    save_results_to_csv(os.path.join(outdir, f"{category}.csv"), result)
    plot_af_with_accuracy(result, os.path.join(outdir, f"{category}"))

    return result


def read_and_process_nipt_samples(chromosome, category, outdir):
    """
    Đọc và xử lý tất cả các mẫu NIPT từ các thư mục trong root_dir.
    """
    total = pd.DataFrame()  # DataFrame tổng

    for ff in PARAMETERS["ff"]:
        result = pd.DataFrame()  # DataFrame cho từng giá trị ff

        for coverage in PARAMETERS["coverage"]:
            for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
                for trio_name, trio_info in TRIO_DATA.items():

                    child = trio_info["child"]
                    mother = trio_info["mother"]
                    father = trio_info["father"]

                    fq_path = os.path.join(
                        fastq_nipt_path(child, mother, father, coverage, ff, index),
                        f"{child}_{mother}_{father}.fastq.gz"
                    )
                    sample_path = os.path.join(
                        statistic_outdir(fq_path, chromosome), 
                        f"{category}.csv"
                    )
                    
                    if not os.path.exists(sample_path):
                        continue

                    # Đọc file CSV
                    df = pd.read_csv(sample_path)

                    # Đảm bảo các cột số được xử lý đúng kiểu
                    numeric_cols = df.columns.difference(['AF (%)'])  # Loại trừ cột 'AF (%)'
                    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

                    # Nếu `result` trống, khởi tạo bằng `df`
                    if result.empty:
                        result = df
                    else:
                        # Cộng dồn dữ liệu
                        result = result.set_index('AF (%)').add(df.set_index('AF (%)'), fill_value=0).reset_index()


        # Thêm cột `ff` vào `result`
        result["ff"] = ff
        result["Child Accuracy"] = df["Glimpse Variants same as Child"] / df["Glimpse Variants"]
        result["Mother Accuracy"] = df["Glimpse Variants same as Mother"] / df["Glimpse Variants"]
        plot_nipt_af_with_accuracy(result, os.path.join(outdir, f"{category}_{ff}.png"))

        # Kết hợp `result` vào `total`
        total = pd.concat([total, result], ignore_index=True)

    save_results_to_csv(os.path.join(outdir, f"{category}_nipt.csv"), total)
    return total

def plot_nipt_af_with_accuracy(df, outdir):
    """
    Vẽ hai biểu đồ (biến thể và độ chính xác) trong cùng một file ảnh.
    
    Parameters:
        df (DataFrame): Dữ liệu từ file CSV.
        outdir (str): Đường dẫn file để lưu biểu đồ.
    """

    # Cột `AF (%)` sẽ là trục X
    x = df["AF (%)"]

    # Các cột cần plot cho biểu đồ biến thể
    columns_to_plot = [
        "Child Variants",
        "Mother Variants",
        "Glimpse Variants",
        "Glimpse Variants same as Child",
        "Glimpse Variants same as Mother",
    ]

    # Tạo figure và chia thành 2 subplot
    plt.figure(figsize=(20, 12))

    # **Biểu đồ 1: Biến thể**
    plt.subplot(2, 1, 1)  # Hàng 1, cột 1 (trên cùng)
    for column in columns_to_plot:
        plt.plot(x, df[column], label=column)

    plt.title("Biểu đồ Biến thể theo AF (%)", fontsize=14)
    plt.xlabel("AF (%)", fontsize=12)
    plt.ylabel("Giá trị", fontsize=12)
    plt.legend(title="Các cột", fontsize=10)
    plt.grid(True)

    # **Biểu đồ 2: Độ chính xác**
    plt.subplot(2, 1, 2)  # Hàng 2, cột 1 (dưới cùng)
    plt.plot(x, df["Child Accuracy"], '--', label="Child Accuracy", color="blue")  # Đường nét đứt màu xanh
    plt.plot(x, df["Mother Accuracy"], '--', label="Mother Accuracy", color="green")  # Đường nét đứt màu xanh lá

    plt.title("Biểu đồ Độ chính xác theo AF (%)", fontsize=14)
    plt.xlabel("AF (%)", fontsize=12)
    plt.ylabel("Độ chính xác", fontsize=12)
    plt.legend(title="Độ chính xác", fontsize=10)
    plt.grid(True)

    # Lưu biểu đồ vào file
    plt.tight_layout()
    plt.savefig(outdir, dpi=300)
    plt.close()

    print(f"Biểu đồ đã được lưu tại: {outdir}")


def plot_af_with_accuracy(df, outdir):
    """
    Vẽ hai biểu đồ (biến thể và độ chính xác) trong cùng một file ảnh.
    
    Parameters:
        df (DataFrame): Dữ liệu từ file CSV.
        outdir (str): Đường dẫn file để lưu biểu đồ.
    """

    # Cột `AF (%)` sẽ là trục X
    x = df["AF (%)"]

    # Các cột cần plot cho biểu đồ biến thể
    columns_to_plot = [
        "Truth Variants",
        "Glimpse Variants",
        "Glimpse True Variants",
    ]

    # Tạo figure và chia thành 2 subplot
    plt.figure(figsize=(20, 12))

    # **Biểu đồ 1: Biến thể**
    plt.subplot(2, 1, 1)  # Hàng 1, cột 1 (trên cùng)
    for column in columns_to_plot:
        plt.plot(x, df[column], label=column)

    plt.title("Biểu đồ Biến thể theo AF (%)", fontsize=14)
    plt.xlabel("AF (%)", fontsize=12)
    plt.ylabel("Giá trị", fontsize=12)
    plt.legend(title="Các cột", fontsize=10)
    plt.grid(True)

    # **Biểu đồ 2: Độ chính xác**
    plt.subplot(2, 1, 2)  # Hàng 2, cột 1 (dưới cùng)
    plt.plot(x, df["Glimpse Accuracy"], '--', label="Glimpse Accuracy", color="orange")  # Đường nét đứt

    plt.title("Biểu đồ Độ chính xác Glimpse theo AF (%)", fontsize=14)
    plt.xlabel("AF (%)", fontsize=12)
    plt.ylabel("Độ chính xác", fontsize=12)
    plt.legend(title="Độ chính xác", fontsize=10)
    plt.grid(True)

    # Lưu biểu đồ vào file
    plt.tight_layout()
    plt.savefig(outdir, dpi=300)
    plt.close()

    print(f"Biểu đồ đã được lưu tại: {outdir}")


if __name__ == "__main__":
    for chromosome in PARAMETERS["chrs"]:
        
        # Đường dẫn để lưu các biểu đồ
        output_dir = os.path.join(PATHS["plot_directory"], chromosome)
        os.makedirs(output_dir, exist_ok=True)

        # Vẽ biểu đồ
        read_and_process_single_samples(chromosome, "summary", output_dir)
        read_and_process_single_samples(chromosome, "rare_summary", output_dir)

        print(f"Plots saved to {output_dir}")

        read_and_process_nipt_samples(chromosome, "summary", output_dir)
        read_and_process_nipt_samples(chromosome, "rare_summary", output_dir)

