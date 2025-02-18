import os, subprocess
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.path_define import fastq_path_lane1, fastq_single_path, fastq_nipt_path
from helper.logger import setup_logger
from helper.file_utils import filter_with_seqtk, filter_and_trim_with_seqtk
from concurrent.futures import ThreadPoolExecutor

logger = setup_logger(os.path.join(PATHS["logs"], "generate.log"))
total_reads = PARAMETERS["refsize"] / PARAMETERS["read_length"]

def generate_random_reads_files(name, coverage, output_prefix):
    input_file = fastq_path_lane1(name)
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    logger.info(f"Generate sample {name} with coverage {coverage}....")

    num_reads = int(coverage * total_reads)

    # Lọc với seqtk
    return filter_and_trim_with_seqtk(input_file, output_file, num_reads)


def generate_merge_files(child_name, mother_name, coverage, ff, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    logger.info(f"Generate nipt sample, child {child_name} and mother {mother_name} with coverage {coverage}....")

    # Tính toán tỷ lệ đọc cho con và mẹ
    child_reads = int(ff * coverage * total_reads)
    mother_reads = int((1 - ff) * coverage * total_reads)

    # Tạo các file FASTQ đã lọc cho con và mẹ
    child_output = f"{output_prefix}_child.fastq.gz"
    mother_output = f"{output_prefix}_mother.fastq.gz"

    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(
            filter_and_trim_with_seqtk,
            [fastq_path_lane1(child_name), fastq_path_lane1(mother_name)],  # Input files
            [child_output, mother_output],  # Output files
            [child_reads, mother_reads]  # num_reads
        )

    # Hợp nhất các tệp con và mẹ
    cmd_merge = f"{TOOLS['zcat']} {child_output} {mother_output} | gzip > {output_file}"
    subprocess.run(cmd_merge, shell=True, check=True)

    try:
        os.remove(child_output)
        os.remove(mother_output)
        logger.info(f"Temporary files {child_output} and {mother_output} deleted.")
    except OSError as e:
        logger.error(f"Error deleting temporary files: {e}")

    return output_file


def generate_single_sample(name, coverage, index):
    """
    Tạo file dữ liệu tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    
    sample_output_dir = fastq_single_path(name, coverage, index)
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, name)
    output_dir = generate_random_reads_files(name, coverage, output_prefix)
    return output_dir


def generate_nipt_sample(child_name, mother_name, father_name, coverage, ff, index):
    """
    Tạo file nipt tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    sample_output_dir = fastq_nipt_path(child_name, mother_name, father_name, coverage, ff, index)
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, f"{child_name}_{mother_name}_{father_name}")
    output_dir = generate_merge_files(child_name, mother_name, coverage, ff, output_prefix)
    return output_dir
 