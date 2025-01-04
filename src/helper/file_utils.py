import gzip
import subprocess
import os
from helper.config import PATHS, TOOLS
from helper.slurm import create_slurm_script, submit_to_slurm

def download_file(url, output_path):
    """
    Hàm tải file từ URL.
    """
    print(f"Downloading file from {url} to {output_path}...")
    result = subprocess.run(["wget", "-O", output_path, url], capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to download file: {result.stderr.decode()}")
    print(f"File downloaded to {output_path}.")

def read_fastq_file(name):
    """
    Đọc file FASTQ theo tên mẫu.
    Lấy đường dẫn từ config PATHS và tên mẫu.
    """
    fastq_path = os.path.join(PATHS["fastq_directory"], f"{name}.fastq.gz")

    print(f"Đọc dữ liệu từ {fastq_path}")

    if not os.path.exists(fastq_path):
        raise FileNotFoundError(f"FASTQ file not found for sample {name} at {fastq_path}")

    headers, reads, plus_separators, qualities = [], [], [], []
    with gzip.open(fastq_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:  # Kết thúc file
                break
            reads.append(f.readline().strip())
            plus_separators.append(f.readline().strip())
            qualities.append(f.readline().strip())
            headers.append(header)

    print(f"Đã đọc dữ liệu FASTQ của mẫu {name} từ {fastq_path}")
    return headers, reads, plus_separators, qualities


def save_to_fastq(output_file, selected_reads, selected_qualities, selected_headers, selected_plus_separators):
    """
    Hàm lưu reads và chất lượng vào file FASTQ.GZ.
    """
    with gzip.open(output_file, 'wt') as f:
        for header, read, plus, quality in zip(selected_headers, selected_reads, selected_plus_separators, selected_qualities):
            f.write(f'{header}\n{read}\n{plus}\n{quality}\n')
    print(f"Đã lưu kết quả vào {output_file}")
    return output_file


def save_results_to_csv(filename, df, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Save the dataframe to a corresponding CSV file
    file_path = os.path.join(output_dir, f"{filename}.csv")
    df.to_csv(file_path, index=False)
    print(f"Saved: {file_path}")


def extract_vcf(sample_name, output_vcf_path, chr=None):
    """
    Tách mẫu VCF từ file tham chiếu bằng cách sử dụng bcftools.
    """
    vcf_reference = PATHS["vcf_reference"]

    # Xây dựng lệnh bcftools
    vcf_command = [
        TOOLS["bcftools"], "view", vcf_reference,
        "--samples", sample_name,
        "-Oz", "-o", output_vcf_path,
    ]

    if chr:
        vcf_command.extend(["--regions", chr])

    # Tạo script SLURM
    script_name = os.path.join(PATHS["vcf_ref"], f"extract_{sample_name}_{chr}.slurm")
    job_name = f"extract_vcf_{sample_name}_{chr}"
    commands = " ".join(vcf_command)

    create_slurm_script(script_name, job_name, commands, PATHS["vcf_directory"])
    job_id = submit_to_slurm(script_name)
    print(f"Submitted job to extract VCF for {sample_name}, chromosome {chr}, Job ID: {job_id}")
    return job_id
