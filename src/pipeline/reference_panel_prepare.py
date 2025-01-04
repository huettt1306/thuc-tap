import os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.slurm import create_slurm_script, submit_to_slurm, wait_for_job_completion
from helper.path_define import vcf_prefix, get_vcf_path, filtered_tsv_path, filtered_vcf_path, chunks_path

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

def download_reference_panel_slurm(chromosome):
    """
    Submit a SLURM job to download reference panel from 1KGP FTP if not already exists.
    """
    url = f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz"
    vcf_path = get_vcf_path(chromosome)
    index_path = f"{vcf_path}.tbi"

    if os.path.exists(vcf_path) and os.path.exists(index_path):
        print(f"=== Reference panel for {chromosome} already exists. Skipping download. ===")
        return vcf_path

    print(f"Preparing SLURM script to download reference panel for chromosome {chromosome}...")
    commands = f"""
wget -c {url} -O {vcf_path}
wget -c {url}.tbi -O {index_path}
"""
    script_name = os.path.join(reference_path, f"download_{chromosome}.slurm")
    job_name = f"download_{chromosome}"

    create_slurm_script(script_name, job_name, commands, reference_path)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted download job for chromosome {chromosome} with Job ID: {job_id}")
    return job_id, vcf_path

def normalize_and_filter_reference_slurm(chromosome):
    """
    Submit a SLURM job to normalize and filter the reference panel.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    output_vcf = os.path.join(reference_path, f"{prefix}.filtered.biallelic.snp.maf0.001.vcf.gz")

    if os.path.exists(output_vcf) and os.path.exists(f"{output_vcf}.tbi"):
        print(f"=== Filtered VCF already exists. Skipping normalization and filtering. ===")
        return output_vcf

    print(f"Preparing SLURM script to normalize and filter VCF {vcf_path}...")
    commands = f"""
{BCFTOOLS} norm -m -any {vcf_path} -Ou --threads 8 | \
    {BCFTOOLS} view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o {output_vcf}
{BCFTOOLS} index -f {output_vcf}
"""
    script_name = os.path.join(reference_path, f"normalize_filter_chr{prefix}.slurm")
    job_name = f"normalize_filter_chr{prefix}"

    create_slurm_script(script_name, job_name, commands, reference_path)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted normalization and filtering job with Job ID: {job_id}")
    return job_id, output_vcf

def process_snp_sites_slurm(chromosome):
    """
    Submit a SLURM job to process SNP sites.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    filtered_vcf = filtered_vcf_path(chromosome)
    tsv_output = filtered_tsv_path(chromosome)

    if os.path.exists(filtered_vcf) and os.path.exists(f"{filtered_vcf}.tbi") and os.path.exists(tsv_output):
        print(f"=== SNP site files already exist. Skipping SNP site processing. ===")
        return filtered_vcf, tsv_output

    print(f"Preparing SLURM script to process SNP sites from VCF {vcf_path}...")
    commands = f"""
{BCFTOOLS} view -G -m 2 -M 2 -v snps {vcf_path} -Oz -o {filtered_vcf} --threads 8
{BCFTOOLS} index -f {filtered_vcf}
{BCFTOOLS} query -f '%CHROM\t%POS\t%REF,%ALT\n' {filtered_vcf} | {BGZIP} -c > {tsv_output}
{TABIX} -s1 -b2 -e2 {tsv_output}
"""
    script_name = os.path.join(reference_path, f"process_sites_chr{prefix}.slurm")
    job_name = f"process_sites_chr{prefix}"

    create_slurm_script(script_name, job_name, commands, reference_path)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted SNP site processing job with Job ID: {job_id}")
    return job_id, filtered_vcf, tsv_output

def chunk_reference_genome_slurm(chromosome, chunk_size=2, buffer_size=2):
    """
    Submit a SLURM job to chunk the reference genome.
    """
    prefix = vcf_prefix(chromosome)
    vcf_path = get_vcf_path(chromosome)
    chunks_output = chunks_path(chromosome)

    if os.path.exists(chunks_output):
        print(f"=== Chunk file already exists: {chunks_output}. Skipping chunking. ===")
        return chunks_output

    print(f"Preparing SLURM script to chunk reference genome for region {chromosome}...")
    commands = f"""
{GLIMPSE_CHUNK} --input {vcf_path} --region {chromosome} \
    --window-mb {chunk_size} --buffer-mb {buffer_size} \
    --output {chunks_output} --sequential
"""
    script_name = os.path.join(reference_path, f"chunk_genome_chr{prefix}.slurm")
    job_name = f"chunk_genome_chr{prefix}"

    create_slurm_script(script_name, job_name, commands, reference_path)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted genome chunking job with Job ID: {job_id}")
    return job_id, chunks_output

def prepare_reference_panel():
    """
    Thực hiện toàn bộ quy trình chuẩn bị reference panel cho các chromosome.
    """
    for chromosome in PARAMETERS["chrs"]:
        os.makedirs(reference_path, exist_ok=True)

        if check_reference_panel(chromosome):
            print(f"=== Reference panel for {chromosome} already exists. Skipping. ===")
            continue

        # Step 1: Download reference panel
        download_job_id, vcf_path = download_reference_panel_slurm(chromosome)
        wait_for_job_completion(download_job_id)

        # Step 2: Normalize and filter reference panel
        normalize_job_id, filtered_vcf = normalize_and_filter_reference_slurm(chromosome)
        wait_for_job_completion(normalize_job_id)

        # Step 3: Process SNP sites
        process_job_id, site_path, tsv_path = process_snp_sites_slurm(chromosome)
        wait_for_job_completion(process_job_id)

        # Step 4: Chunk reference genome
        chunk_job_id, chunk_file = chunk_reference_genome_slurm(chromosome)
        wait_for_job_completion(chunk_job_id)

        print(f"=== Reference panel preparation completed for {chromosome} ===")
