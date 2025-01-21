import os, subprocess
from helper.config import PATHS, TOOLS
from helper.path_define import fastq_path, fastq_single_path, fastq_nipt_path, cram_path
from helper.logger import setup_logger
from helper.converter import convert_cram_to_fastq
from concurrent.futures import ThreadPoolExecutor

logger = setup_logger(os.path.join(PATHS["logs"], "generate.log"))

generated_files = []

def filter_with_seqtk(name, output_file, fraction):
    """
    Sử dụng seqtk để lấy mẫu dữ liệu từ file FASTQ theo tỷ lệ nhất định.
    """
    logger.info(f"Filtering with seqtk sample {name} with fraction {fraction}....")
    input_file = fastq_path(name)
    if not os.path.exists(input_file):
        convert_cram_to_fastq(cram_path(name), input_file) 
        if not os.path.exists(input_file):
            logger.error(f"Sample {name} cannot read.")
            raise RuntimeError(f"Failed to read sample: {name}")

    cmd = f"{TOOLS['seqtk']} sample {input_file} {fraction} | gzip > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Filter {name} done.")
    return output_file

def generate_random_reads_files(name, coverage, avg_coverage, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    logger.info(f"Generate sample {name} with coverage {coverage}....")

    # Tính tỷ lệ đọc cần trích xuất
    ratio = coverage / avg_coverage

    # Lọc với seqtk
    return filter_with_seqtk(name, output_file, ratio)


def generate_merge_files(child_name, mother_name, coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        logger.info(f"File {output_file} already exists. Skipping creation.")
        return output_file

    logger.info(f"Generate nipt sample, child {child_name} and mother {mother_name} with coverage {coverage}....")

    # Tính toán tỷ lệ đọc cho con và mẹ
    child_ratio = ff * coverage / child_avg_coverage
    mother_ratio = (1 - ff) * coverage / mother_avg_coverage

    # Tạo các file FASTQ đã lọc cho con và mẹ
    child_output = f"{output_prefix}_child.fastq.gz"
    mother_output = f"{output_prefix}_mother.fastq.gz"

    with ThreadPoolExecutor(max_workers=2) as executor:
        # Submit các tác vụ lọc cho con và mẹ
        future_child = executor.submit(filter_with_seqtk, child_name, child_output, child_ratio)
        future_mother = executor.submit(filter_with_seqtk, mother_name, mother_output, mother_ratio)

        # Đợi các tác vụ hoàn thành
        future_child.result()
        future_mother.result()

    # Hợp nhất các tệp con và mẹ mà không xáo trộn
    cmd_merge = f"{TOOLS['zcat']} {child_output} {mother_output} | gzip > {output_file}"
    subprocess.run(cmd_merge, shell=True, check=True)

    try:
        os.remove(child_output)
        os.remove(mother_output)
        logger.info(f"Temporary files {child_output} and {mother_output} deleted.")
    except OSError as e:
        logger.error(f"Error deleting temporary files: {e}")

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
    return output_dir


def generate_nipt_sample(child_name, mother_name, father_name, child_avg_coverage, mother_avg_coverage, coverage, ff, index):
    """
    Tạo file nipt tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    sample_output_dir = fastq_nipt_path(child_name, mother_name, father_name, coverage, ff, index)
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, f"{child_name}_{mother_name}_{father_name}")
    output_dir = generate_merge_files(child_name, mother_name, coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix)
    return output_dir
 