import gzip
import subprocess
import os
import pandas as pd
from cyvcf2 import VCF
from helper.config import PATHS, TOOLS, PARAMETERS
from helper.logger import setup_logger
from helper.variant import convert_genotype

logger = setup_logger(os.path.join(PATHS["logs"], "file_utils.log"))

def extract_land1_fq(fq_path, r1_path):
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
    
    cmd = f"{TOOLS['seqtk']} seq -1 {fq_path} | gzip > {r1_path}"
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Extracted land 1 for {fq_path}")


def filter_with_seqtk(input_file, output_file, fraction):
    """
    Sử dụng seqtk để lấy mẫu dữ liệu từ file FASTQ theo tỷ lệ nhất định.
    """
    logger.info(f"Filtering with seqtk sample {input_file} with fraction {fraction}....")
    if not os.path.exists(input_file):
        logger.error(f"Sample {input_file} cannot read.")
        raise RuntimeError(f"Failed to read sample: {input_file}")

    cmd = f"{TOOLS['seqtk']} sample {input_file} {fraction} | gzip > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Filter {input_file} done.")
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