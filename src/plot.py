import os
import pandas as pd
import matplotlib.pyplot as plt
from helper.config import PATHS, PARAMETERS, TRIO_DATA
from helper.path_define import statistic_outdir, statistic_nipt_outdir, fastq_single_path, fastq_nipt_path

def num_variants(df, category, method):
    return df.loc[(df['Category'] == category) & (df['Method'] == method), f'{method} Variants'].values[0]

def num_variants_true(df, category, method):
    return df.loc[(df['Category'] == category) & (df['Method'] == method), "True Positives"].values[0]

# Hàm đọc và xử lý dữ liệu
def read_and_process_single_samples(chromosome, category):
    """
    Đọc và xử lý tất cả các mẫu từ các thư mục trong root_dir.
    """
    stats_list = []
    for coverage in PARAMETERS["coverage"]:

        basevar_variants = 0
        glimpse_variants = 0
        basevar_true = 0
        glimpse_true = 0

        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            for trio_name, trio_info in TRIO_DATA.items():
                for role, name in trio_info.items():
                    fq_path = os.path.join(fastq_single_path(name, coverage, index), f"{name}.fastq.gz")
                    sample_path = os.path.join(statistic_outdir(fq_path, chromosome), "summary.csv")                   
                    if not os.path.exists(sample_path):
                        continue

                    df = pd.read_csv(sample_path)
                    
                    basevar_variants += num_variants(df, category, "BaseVar")
                    glimpse_variants += num_variants(df, category, "Glimpse")
                    basevar_true += num_variants_true(df, category, "BaseVar")
                    glimpse_true += num_variants_true(df, category, "Glimpse")

        basevar_ratio = basevar_true / basevar_variants
        glimpse_ratio = glimpse_true / glimpse_variants


        stats_list.append({
            "Coverage": coverage,
            "BaseVar_Ratio": basevar_ratio,
            "Glimpse_Ratio": glimpse_ratio,
        })
    
    return pd.DataFrame(stats_list)

def read_and_process_nipt_samples(chromosome, category):
    """
    Đọc và xử lý tất cả các mẫu NIPT từ các thư mục trong root_dir.
    """
    stats_list = []
    for ff in PARAMETERS["ff"]:
        # Khởi tạo các biến để lưu số liệu
        basevar_variants = {"child": 0, "mother": 0}
        glimpse_variants = {"child": 0, "mother": 0}
        basevar_true = {"child": 0, "mother": 0}
        glimpse_true = {"child": 0, "mother": 0}

        for coverage in PARAMETERS["coverage"]:
            for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
                for trio_name, trio_info in TRIO_DATA.items():

                    child = trio_info["child"]
                    mother = trio_info["mother"]
                    father = trio_info["father"]

                    # Đường dẫn đến dữ liệu của con và mẹ
                    sample = os.path.join(fastq_nipt_path(child, mother, father, coverage, ff, index), f"{child}_{mother}_{father}.fastq.gz")
                    child_path = os.path.join(statistic_nipt_outdir(sample, chromosome, "child"), "summary.csv")
                    mother_path = os.path.join(statistic_nipt_outdir(sample, chromosome, "mother"), "summary.csv")


                    # Đọc dữ liệu con
                    if os.path.exists(child_path):
                        child_df = pd.read_csv(child_path)
                        basevar_variants["child"] += num_variants(child_df, category, "BaseVar")
                        glimpse_variants["child"] += num_variants(child_df, category, "Glimpse")
                        basevar_true["child"] += num_variants_true(child_df, category, "BaseVar")
                        glimpse_true["child"] += num_variants_true(child_df, category, "Glimpse")

                    # Đọc dữ liệu mẹ
                    if os.path.exists(mother_path):
                        mother_df = pd.read_csv(mother_path)
                        basevar_variants["mother"] += num_variants(mother_df, category, "BaseVar")
                        glimpse_variants["mother"] += num_variants(mother_df, category, "Glimpse")
                        basevar_true["mother"] += num_variants_true(mother_df, category, "BaseVar")
                        glimpse_true["mother"] += num_variants_true(mother_df, category, "Glimpse")

        # Tính tỷ lệ
        basevar_ratio_child = basevar_true["child"] / basevar_variants["child"] if basevar_variants["child"] > 0 else 0
        glimpse_ratio_child = glimpse_true["child"] / glimpse_variants["child"] if glimpse_variants["child"] > 0 else 0
        basevar_ratio_mother = basevar_true["mother"] / basevar_variants["mother"] if basevar_variants["mother"] > 0 else 0
        glimpse_ratio_mother = glimpse_true["mother"] / glimpse_variants["mother"] if glimpse_variants["mother"] > 0 else 0

        print(basevar_variants["child"])
        print(basevar_variants["mother"])

        # Thêm vào danh sách kết quả
        stats_list.append({
            "FF": ff,
            "Coverage": coverage,
            "Child_BaseVar_Ratio": basevar_ratio_child,
            "Child_Glimpse_Ratio": glimpse_ratio_child,
            "Mother_BaseVar_Ratio": basevar_ratio_mother,
            "Mother_Glimpse_Ratio": glimpse_ratio_mother,
        })

    return pd.DataFrame(stats_list)



# Hàm vẽ biểu đồ
def plot_data(df, category, output_dir):
    """
    Vẽ biểu đồ tương quan giữa coverage và ratio.
    """
    output_file = os.path.join(output_dir, f"{category}_Coverage_vs_Ratio.png")
    plt.figure(figsize=(12, 6))

    # Biểu đồ cho Total
    plt.subplot(1, 2, 1)
    plt.plot(df['Coverage'], df['BaseVar_Ratio'], label=f"{category} BaseVar", marker='o', color='blue')
    plt.plot(df['Coverage'], df['Glimpse_Ratio'], label=f"{category} Glimpse", marker='o', color='red')
    plt.title(f"{category} Coverage vs Ratio")
    plt.xlabel("Coverage")
    plt.ylabel("Ratio")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Biểu đồ đã được lưu tại: {output_file}")


def plot_data_by_ff(df, category, output_dir):
    """
    Vẽ biểu đồ tương quan giữa ff và ratio cho các giá trị tính được.
    """
    output_file = os.path.join(output_dir, f"{category}_FF_vs_Ratio.png")
    plt.figure(figsize=(12, 6))

    # Biểu đồ cho Child
    plt.subplot(1, 2, 1)
    plt.plot(df['FF'], df['Child_BaseVar_Ratio'], label=f"Child BaseVar", marker='o', color='blue')
    plt.plot(df['FF'], df['Child_Glimpse_Ratio'], label=f"Child Glimpse", marker='o', color='green')
    plt.title(f"{category} FF vs Ratio (Child)")
    plt.xlabel("FF")
    plt.ylabel("Ratio")
    plt.legend()
    plt.grid()

    # Biểu đồ cho Mother
    plt.subplot(1, 2, 2)
    plt.plot(df['FF'], df['Mother_BaseVar_Ratio'], label=f"Mother BaseVar", marker='o', color='red')
    plt.plot(df['FF'], df['Mother_Glimpse_Ratio'], label=f"Mother Glimpse", marker='o', color='purple')
    plt.title(f"{category} FF vs Ratio (Mother)")
    plt.xlabel("FF")
    plt.ylabel("Ratio")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Biểu đồ đã được lưu tại: {output_file}")


if __name__ == "__main__":
    for chromosome in PARAMETERS["chrs"]:
        
        # Đường dẫn để lưu các biểu đồ
        output_dir = os.path.join(PATHS["plot_directory"], chromosome)
        os.makedirs(output_dir, exist_ok=True)

        # Vẽ biểu đồ
        plot_data(read_and_process_single_samples(chromosome, "Total"), "Total", output_dir)
        plot_data(read_and_process_single_samples(chromosome, "Rare"), "Rare", output_dir)

        print(f"Plots saved to {output_dir}")

        plot_data_by_ff(read_and_process_nipt_samples(chromosome, "Total"), "Total", output_dir)
        plot_data_by_ff(read_and_process_nipt_samples(chromosome, "Rare"), "Rare", output_dir)

