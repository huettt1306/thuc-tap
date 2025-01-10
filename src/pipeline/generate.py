import os, subprocess
import random
from helper.file_utils import save_to_fastq, read_fastq_file
from helper.config import PATHS, TOOLS
from helper.path_define import fastq_path, fastq_single_path, fastq_nipt_path
from helper.logger import setup_logger

logger = setup_logger(os.path.join(PATHS["logs"], "generate.log"))

generated_files = []

def filter_with_seqtk(input_file, output_file, fraction):
    """
    Sử dụng seqtk để lấy mẫu dữ liệu từ file FASTQ theo tỷ lệ nhất định.
    """
    cmd = f"{TOOLS['seqtk']} sample -s100 {input_file} {fraction} | gzip > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    return output_file

def generate_random_reads_files(name, coverage, avg_coverage, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    # Tính tỷ lệ đọc cần trích xuất
    ratio = coverage / avg_coverage

    # Lọc với seqtk
    return filter_with_seqtk(fastq_path(name), output_file, ratio)


def generate_merge_files(child_name, mother_name, coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    # Tính toán tỷ lệ đọc cho con và mẹ
    child_ratio = ff * coverage / child_avg_coverage
    mother_ratio = (1 - ff) * coverage / mother_avg_coverage

    # Tạo các file FASTQ đã lọc cho con và mẹ
    child_output = f"{output_prefix}_child.fastq.gz"
    mother_output = f"{output_prefix}_mother.fastq.gz"

    filter_with_seqtk(fastq_path(child_name), child_output, child_ratio)
    filter_with_seqtk(fastq_path(mother_name), mother_output, mother_ratio)

    # Hợp nhất các tệp con và mẹ
    cmd_merge = f"{TOOLS['zcat']} {child_output} {mother_output} | shuf | gzip > {output_file}"
    subprocess.run(cmd_merge, shell=True, check=True)

    return output_file


def generate_single_sample(name, avg_coverage, coverage, index):
    """
    Tạo file dữ liệu tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    
    sample_output_dir = fastq_single_path(name, coverage, index)
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, name)
    output_dir = generate_random_reads_files(name, coverage, avg_coverage, output_prefix)
    generated_files.append(output_dir)
    return output_dir


def generate_nipt_sample(child_name, child_avg_coverage, mother_name, mother_avg_coverage, coverage, ff, index):
    """
    Tạo file nipt tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    sample_output_dir = fastq_nipt_path(child_name, mother_name, coverage, ff, index)
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, f"{child_name}_{mother_name}")
    output_dir = generate_merge_files(child_name, mother_name, coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix)
    generated_files.append(output_dir)
    return output_dir
 