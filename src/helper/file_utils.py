import gzip
import subprocess
import os
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.logger import setup_logger

logger = setup_logger(os.path.join(PATHS["logs"], "file_utils.log"))


def download_file(url, output_path):
    """
    Hàm tải file từ URL.
    """
    logger.info(f"Downloading file from {url} to {output_path} ")
    command = ["wget", "-O", output_path, url]

    try:
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to download file: {result.stderr}")
            raise RuntimeError(f"Failed to download file: {result.stderr}")
        logger.info(f"File downloaded successfully to {output_path}.")
    except Exception as e:
        logger.error(f"Error during file download: {e}")
        raise

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
            sequence = f.readline().strip()
            plus_separator = f.readline().strip()
            quality = f.readline().strip()

            # Kiểm tra số dòng đọc được có hợp lệ không (4 dòng cho mỗi read)
            if not header or not sequence or not plus_separator or not quality:
                logger.error(f"File {fastq_path} không hợp lệ: thiếu dữ liệu cho một hoặc nhiều reads.")
                break

            headers.append(header)
            reads.append(sequence)
            plus_separators.append(plus_separator)
            qualities.append(quality)

    print(f"Đã đọc dữ liệu FASTQ của mẫu {name} từ {fastq_path}")
    return reads, qualities, headers, plus_separators


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
    threads = PARAMETERS["bcftools"]["threads"]
    vcf_reference = PATHS["vcf_reference"]

    # Xây dựng lệnh bcftools
    vcf_command = [
        TOOLS["bcftools"], "view", vcf_reference,
        "--samples", sample_name,
        "-Oz", "-o", output_vcf_path,
        f"--threads={threads}"
    ]

    if chr:
        vcf_command.extend(["--regions", chr])

    try:
        result = subprocess.run(vcf_command, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to extract VCF: {result.stderr}")
            raise RuntimeError(f"Failed to extract VCF: {result.stderr}")
        logger.info(f"VCF file extracted successfully to {output_vcf_path}.")
    except Exception as e:
        logger.error(f"Error extracting VCF: {e}")
        raise

