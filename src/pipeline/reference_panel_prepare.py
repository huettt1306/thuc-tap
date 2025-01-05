import subprocess
import os
import logging
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import vcf_prefix, get_vcf_path, filtered_tsv_path, filtered_vcf_path, chunks_path
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(log_file="logs/reference_panel_pipeline.log")

BCFTOOLS = TOOLS["bcftools"]
BGZIP = TOOLS["bgzip"]
TABIX = TOOLS["tabix"]
GLIMPSE_CHUNK = TOOLS["GLIMPSE_chunk"]
reference_path = PATHS["reference_path"]

def check_reference_panel(chromosome):
    """
    Kiểm tra dữ liệu reference panel đã tồn tại hay chưa.
    """
    prefix = vcf_prefix(chromosome)
    required_files = [
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.vcf.gz",
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.vcf.gz.tbi",
        f"{reference_path}/{prefix}.biallelic.snp.maf0.001.sites.tsv.gz",
        f"{reference_path}/{prefix}.chunks.txt"
    ]
    for file in required_files:
        if not os.path.exists(file):
            return False
    return True

def download_reference_panel(chromosome):
    """
    Download reference panel from 1KGP FTP if not already exists.
    """
    url = f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz"
    vcf_path = get_vcf_path(chromosome)
    index_path = f"{vcf_path}.tbi"

    if os.path.exists(vcf_path) and os.path.exists(index_path):
        logger.info(f"Reference panel for {chromosome} already exists. Skipping download.")
        return vcf_path

    logger.info(f"Downloading reference panel for chromosome {chromosome}...")
    commands = [
        ["wget", "-c", url, "-O", vcf_path],
        ["wget", "-c", f"{url}.tbi", "-O", index_path]
    ]
    
    for command in commands:
        process = subprocess.run(command, capture_output=True, text=True)
        if process.returncode != 0:
            logger.error(f"Error downloading file: {process.stderr}")
            raise RuntimeError(f"Error downloading file: {process.stderr}")

    logger.info(f"Downloaded reference panel for chromosome {chromosome}.")
    return vcf_path

def normalize_and_filter_reference(chromosome):
    """
    Normalize and filter the reference panel.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    output_vcf = os.path.join(reference_path, f"{prefix}.filtered.biallelic.snp.maf0.001.vcf.gz")

    if os.path.exists(output_vcf) and os.path.exists(f"{output_vcf}.tbi"):
        logger.info(f"Filtered VCF already exists. Skipping normalization and filtering.")
        return output_vcf

    logger.info(f"Normalizing and filtering VCF for chromosome {chromosome}...")
    command = [
        BCFTOOLS, "norm", "-m", "-any", vcf_path, "-Ou",
        "--threads", "1", "|",
        BCFTOOLS, "view", "-m", "2", "-M", "2", "-v", "snps", "-i", "'MAF>0.001'",
        "--threads", "1", "-Oz", "-o", output_vcf
    ]

    process = subprocess.run(" ".join(command), shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error normalizing and filtering: {process.stderr}")
        raise RuntimeError(f"Error normalizing and filtering: {process.stderr}")

    index_command = [TABIX, "-f", output_vcf]
    subprocess.run(index_command, check=True)

    logger.info(f"Filtered VCF created at {output_vcf}.")
    return output_vcf

def process_snp_sites(chromosome):
    """
    Process SNP sites.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    filtered_vcf = filtered_vcf_path(chromosome)
    tsv_output = filtered_tsv_path(chromosome)

    if os.path.exists(filtered_vcf) and os.path.exists(f"{filtered_vcf}.tbi") and os.path.exists(tsv_output):
        logger.info(f"SNP site files already exist. Skipping SNP site processing.")
        return filtered_vcf, tsv_output

    logger.info(f"Processing SNP sites for chromosome {chromosome}...")
    commands = [
        [BCFTOOLS, "view", "-G", "-m", "2", "-M", "2", "-v", "snps", vcf_path, "-Oz", "-o", filtered_vcf, "--threads", "1"],
        [TABIX, "-f", filtered_vcf],
        [BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF,%ALT\n", filtered_vcf, "|", BGZIP, "-c", "-o", tsv_output],
        [TABIX, "-s1", "-b2", "-e2", tsv_output]
    ]

    for command in commands:
        process = subprocess.run(command, capture_output=True, text=True, shell=True)
        if process.returncode != 0:
            logger.error(f"Error processing SNP sites: {process.stderr}")
            raise RuntimeError(f"Error processing SNP sites: {process.stderr}")

    logger.info(f"Processed SNP sites for chromosome {chromosome}.")
    return filtered_vcf, tsv_output

def chunk_reference_genome(chromosome, chunk_size=2, buffer_size=2):
    """
    Chunk the reference genome.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    chunks_output = chunks_path(chromosome)

    if os.path.exists(chunks_output):
        logger.info(f"Chunk file already exists: {chunks_output}. Skipping chunking.")
        return chunks_output

    logger.info(f"Chunking reference genome for chromosome {chromosome}...")
    command = [
        GLIMPSE_CHUNK, "--input", vcf_path, "--region", chromosome,
        "--window-mb", str(chunk_size), "--buffer-mb", str(buffer_size),
        "--output", chunks_output, "--sequential"
    ]

    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error chunking reference genome: {process.stderr}")
        raise RuntimeError(f"Error chunking reference genome: {process.stderr}")

    logger.info(f"Chunk file created at {chunks_output}.")
    return chunks_output

def prepare_reference_panel():
    """
    Thực hiện toàn bộ quy trình chuẩn bị reference panel cho các chromosome.
    """
    for chromosome in PARAMETERS["chrs"]:
        os.makedirs(reference_path, exist_ok=True)

        if check_reference_panel(chromosome):
            logger.info(f"Reference panel for {chromosome} already exists. Skipping.")
            continue

        # Step 1: Download reference panel
        download_reference_panel(chromosome)

        # Step 2: Normalize and filter reference panel
        normalize_and_filter_reference(chromosome)

        # Step 3: Process SNP sites
        process_snp_sites(chromosome)

        # Step 4: Chunk reference genome
        chunk_reference_genome(chromosome)

        logger.info(f"Reference panel preparation completed for {chromosome}.")
