import gzip
import subprocess
import os
import pandas as pd
from cyvcf2 import VCF
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.logger import setup_logger
from helper.converter import convert_cram_to_fastq
from helper.variant import convert_genotype

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

    logger.info(f"Đọc dữ liệu từ {fastq_path}")

    if not os.path.exists(fastq_path):
        cram_path = os.path.join(PATHS["cram_directory"], f"{name}.final.cram")
        convert_cram_to_fastq(cram_path, fastq_path) 
        if not os.path.exists(fastq_path):
            logger.error(f"Sample {name} cannot read.")
            raise RuntimeError(f"Failed to read sample: {name}")


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

    logger.info(f"Đã đọc dữ liệu FASTQ của mẫu {name} từ {fastq_path}")
    return reads, qualities, headers, plus_separators


def save_to_fastq(output_file, selected_reads, selected_qualities, selected_headers, selected_plus_separators):
    """
    Hàm lưu reads và chất lượng vào file FASTQ.GZ.
    """
    with gzip.open(output_file, 'wt') as f:
        for header, read, plus, quality in zip(selected_headers, selected_reads, selected_plus_separators, selected_qualities):
            f.write(f'{header}\n{read}\n{plus}\n{quality}\n')
    logger.info(f"Đã lưu kết quả vào {output_file}")
    return output_file


def save_results_to_csv(file_path, df):
    # Save the dataframe to the specified file path
    os.makedirs(os.path.dirname(file_path), exist_ok=True)  # Tạo thư mục nếu chưa tồn tại

    df.to_csv(file_path, index=False)
    logger.info(f"Saved: {file_path}")


def extract_vcf(sample_name, vcf_reference, output_vcf_path):
    """
    Tách mẫu VCF từ file tham chiếu bằng cách sử dụng bcftools.
    """

    # Xây dựng lệnh bcftools
    vcf_command = [
        TOOLS["bcftools"], "view", vcf_reference,
        "--samples", sample_name,
        "-Oz", "-o", output_vcf_path,
        f"--threads={PARAMETERS['threads']}"
    ]

    try:
        result = subprocess.run(vcf_command, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to extract VCF: {result.stderr}")
            raise RuntimeError(f"Failed to extract VCF: {result.stderr}")
        logger.info(f"VCF file extracted successfully to {output_vcf_path}.")
    except Exception as e:
        logger.error(f"Error extracting VCF: {e}")
        raise

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