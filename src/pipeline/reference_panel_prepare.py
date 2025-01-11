import subprocess
import os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import vcf_prefix, get_vcf_path, filtered_tsv_path, filtered_vcf_path, chunks_path, norm_vcf_path
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "reference_panel_pipeline.log"))


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
        ["wget", "-c", url, "-P", reference_path],
        ["wget", "-c", f"{url}.tbi", "-P", reference_path]
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
    vcf_path = get_vcf_path(chromosome)
    output_vcf = norm_vcf_path(chromosome)

    if os.path.exists(output_vcf) and os.path.exists(f"{output_vcf}.tbi"):
        logger.info(f"Filtered VCF already exists. Skipping normalization and filtering.")
        return output_vcf

    logger.info(f"Normalizing and filtering VCF for chromosome {chromosome}...")
    command = [
        BCFTOOLS, "norm", "-m", "-any", vcf_path, "-Ou",
        "--threads", f"{PARAMETERS['threads']}", "|",
        BCFTOOLS, "view", "-m", "2", "-M", "2", "-v", "snps", "-i", "'MAF>0.001'",
        "--threads", f"{PARAMETERS['threads']}", "-Oz", "-o", output_vcf
    ]

    process = subprocess.run(" ".join(command), shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error normalizing and filtering: {process.stderr}")
        raise RuntimeError(f"Error normalizing and filtering: {process.stderr}")

    index_command = [BCFTOOLS, "index", "-f", output_vcf]
    subprocess.run(index_command, check=True)

    logger.info(f"Filtered VCF created at {output_vcf}.")
    return output_vcf

def process_snp_sites(chromosome):
    """
    Process SNP sites.
    """
    vcf_path = norm_vcf_path(chromosome)
    filtered_vcf = filtered_vcf_path(chromosome)
    tsv_output = filtered_tsv_path(chromosome)

    if os.path.exists(filtered_vcf) and os.path.exists(f"{filtered_vcf}.tbi") and os.path.exists(tsv_output):
        logger.info(f"SNP site files already exist. Skipping SNP site processing.")
        return filtered_vcf, tsv_output

    logger.info(f"Processing SNP sites for chromosome {chromosome}...")
    commands = [
        [BCFTOOLS, "view", "-G", "-m", "2", "-M", "2", "-v", "snps", vcf_path, "-Oz", "-o", filtered_vcf, "--threads", f"{PARAMETERS['threads']}",],
        [BCFTOOLS, "index", "-f", "-@", f"{PARAMETERS['threads']}", filtered_vcf],
    ]

    # Chạy các lệnh bcftools view và index
    for command in commands:
        process = subprocess.run(command, check=True)
        if process.returncode != 0:
            logger.error(f"Error processing SNP sites: {process.stderr}")
            raise RuntimeError(f"Error processing SNP sites: {process.stderr}")

    # Chạy bcftools query và bgzip
    logger.info(f"Running bcftools query and bgzip for chromosome {chromosome}...")
    with open(tsv_output, "wb") as output_file:
        query_command = [BCFTOOLS, "query", "-f", "%CHROM\\t%POS\\t%REF,%ALT\\n", filtered_vcf]
        bgzip_command = [BGZIP, "-c"]
        query_process = subprocess.Popen(query_command, stdout=subprocess.PIPE)
        subprocess.run(bgzip_command, stdin=query_process.stdout, stdout=output_file, check=True)
        query_process.stdout.close()
        query_process.wait()

    # Chạy tabix để tạo index cho tsv_output
    logger.info(f"Running tabix for chromosome {chromosome}...")
    tabix_command = [TABIX, "-s1", "-b2", "-e2", "-@", f"{PARAMETERS['threads']}", tsv_output]
    subprocess.run(tabix_command, check=True)

    logger.info(f"Processed SNP sites for chromosome {chromosome}.")
    return filtered_vcf, tsv_output


def chunk_reference_genome(chromosome):
    """
    Chunk the reference genome.
    """
    vcf_path = get_vcf_path(chromosome)
    chunks_output = chunks_path(chromosome)

    if os.path.exists(chunks_output):
        logger.info(f"Chunk file already exists: {chunks_output}. Skipping chunking.")
        return chunks_output

    logger.info(f"Chunking reference genome for chromosome {chromosome}...")
    command = [
        GLIMPSE_CHUNK, "--input", vcf_path, "--region", chromosome,
        "--window-mb", "2", "--buffer-mb", "0.2",
        "--output", chunks_output, "--sequential"
    ]

    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error chunking reference genome: {process.stderr}")
        raise RuntimeError(f"Error chunking reference genome: {process.stderr}")

    logger.info(f"Chunk file created at {chunks_output}.")
    return chunks_output

def prepare_gatk_bundle():
    dbsnp = os.path.join(PATHS["gatk_bundle_dir"], "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),

    if not os.path.exists(dbsnp):
        print(f"{dbsnp} not found. Compressing {dbsnp}...")
        
        # Sử dụng bgzip để nén tệp .vcf thành .vcf.gz
        bgzip_cmd = [TOOLS['bgzip'], dbsnp]
        subprocess.run(bgzip_cmd, check=True)
        
        # Lập chỉ mục tệp .vcf.gz bằng tabix
        tabix_cmd = [TOOLS['tabix'], "-f", "-@", f"{PARAMETERS['threads']}", f"{dbsnp}.gz"]
        subprocess.run(tabix_cmd, check=True)
        
        print(f"{dbsnp} has been compressed and indexed.")
    else:
        print(f"{dbsnp} already exists. No action needed.")


def prepare_reference_panel():
    """
    Thực hiện toàn bộ quy trình chuẩn bị reference panel cho các chromosome.
    """
    for chromosome in PARAMETERS["chrs"]:
        os.makedirs(reference_path, exist_ok=True)

        if check_reference_panel(chromosome):
            logger.info(f"Reference panel for {chromosome} already exists. Skipping.")
            continue
        
        # Step 0: verify gatk bundle
        prepare_gatk_bundle()

        # Step 1: Download reference panel
        download_reference_panel(chromosome)

        # Step 2: Normalize and filter reference panel
        normalize_and_filter_reference(chromosome)

        # Step 3: Process SNP sites
        process_snp_sites(chromosome)

        # Step 4: Chunk reference genome
        chunk_reference_genome(chromosome)

        logger.info(f"Reference panel preparation completed for {chromosome}.")
