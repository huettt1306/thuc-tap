import random
import subprocess
import os
import pandas as pd
from cyvcf2 import VCF
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.logger import setup_logger
from helper.converter import convert_genotype
from statistic.ALT import valid_alt

logger = setup_logger(os.path.join(PATHS["logs"], "file_utils.log"))

def extract_lane1_fq(fq_path, r1_path):
    """
    Sử dụng seqtk để lấy mẫu dữ liệu từ file FASTQ theo tỷ lệ nhất định.
    """
    logger.info(f"Extracting land 1 for {fq_path}...")
    if os.path.exists(r1_path):
        logger.info("Land 1 fastq already exists.")
        return
    if not os.path.exists(fq_path):
        logger.error(f"Sample {fq_path} cannot read.")
        raise RuntimeError(f"Failed to read sample: {fq_path}")
    
    cmd = f'{TOOLS["seqkit"]} grep -rp ":1:" {fq_path} -o {r1_path}'
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Extracted land 1 for {fq_path}")


def filter_with_seqtk(input_file, output_file, fraction):
    """
    Sử dụng seqtk để lấy mẫu dữ liệu từ file FASTQ theo tỷ lệ nhất định.
    """
    logger.info(f"Filtering with seqtk sample {input_file} with fraction {fraction}....")
    print(f"Filtering with seqtk sample {input_file} with fraction {fraction}....")
    if not os.path.exists(input_file):
        logger.error(f"Sample {input_file} cannot be read.")
        raise RuntimeError(f"Failed to read sample: {input_file}")

    seed = seed = random.randint(1, 10**9 + 7)

    cmd = f"{TOOLS['seqtk']} sample -s {seed} {input_file} {fraction} | gzip > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Filter {input_file} done.")
    return output_file

def filter_and_trim_with_seqtk(input_file, output_file, num_reads, max_length=PARAMETERS["read_length"]):
    logger.info(f"Filtering {num_reads} reads from {input_file} and trimming to {max_length}bp...")
    print(f"Filtering {num_reads} reads from {input_file} and trimming to {max_length}bp...")

    temp_output = os.path.join(os.path.dirname(output_file), "tmp.fastq")

    if not os.path.exists(input_file):
        logger.error(f"Sample {input_file} cannot be read.")
        raise RuntimeError(f"Failed to read sample: {input_file}")

    seed = random.randint(1, 10**9 + 7)

    cmd_sample = f"{TOOLS['seqtk']} sample -s {seed} {input_file} {num_reads} > {temp_output}"
    subprocess.run(cmd_sample, shell=True, check=True)

    cmd_trim = f"{TOOLS['seqtk']} trimfq -L {max_length} {temp_output} | gzip > {output_file}"
    subprocess.run(cmd_trim, shell=True, check=True)

    os.remove(temp_output)

    logger.info(f"Filtering and trimming done. Output saved to {output_file}.")
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
        "-m", "2", "-M", "2",
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

def process_vcf(vcf_path, method_name="Test"):
    """
    Đọc VCF và trích xuất thông tin cần thiết.
    """
    variants = []
    try:
        logger.info(f"Processing VCF file: {vcf_path} for method {method_name}")
        vcf_reader = VCF(vcf_path)
        for record in vcf_reader:
            variants.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT": record.ALT[0],
                f"AF_{method_name}": dict(record.INFO).get('AF', -1),
                f"GT_{method_name}": convert_genotype(record.genotypes[0]),
                method_name: True
            })

    except Exception as e:
        logger.error(f"Error processing VCF file {vcf_path}: {e}")
        raise

    logger.info(f"Finished processing VCF file: {vcf_path}")
    return pd.DataFrame(variants)