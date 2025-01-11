import subprocess
import os, re
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import bamlist_dir, glimpse_outdir
from helper.path_define import filtered_vcf_path, filtered_tsv_path, chunks_path, norm_vcf_path
from helper.logger import setup_logger

# Thiết lập logger
logger = setup_logger(os.path.join(PATHS["logs"], "glimpse_pipeline.log"))


BCFTOOLS = TOOLS["bcftools"]
BGZIP = TOOLS["bgzip"]
TABIX = TOOLS["tabix"]
GLIMPSE_PHASE = TOOLS["GLIMPSE_phase"]
GLIMPSE_LIGATE = TOOLS["GLIMPSE_ligate"]
REF = PATHS["ref"]
MAP_PATH = PATHS["map_path"]

def compute_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    os.makedirs(glpath, exist_ok=True)

    bamlist = bamlist_dir(fq)
    with open(bamlist, "r") as bam_file:
        bam_lines = bam_file.readlines()

    for line in bam_lines:
        line = line.strip()
        filename = os.path.basename(line)
        name = filename.split(".")[0]

        output_vcf = os.path.join(glpath, f"{name}.{chromosome}.vcf.gz")
        command = [
            BCFTOOLS, "mpileup", "-@", f"{PARAMETERS['threads']}",
            "-f", REF, "-I", "-E", "-a", "FORMAT/DP",
            "-T", filtered_vcf_path(chromosome), "-r", chromosome, line, "-Ou"
        ]
        call_command = [
            BCFTOOLS, "call", "-Aim", "-C", "alleles", "-@", f"{PARAMETERS['threads']}",
            "-T", filtered_tsv_path(chromosome), "-Oz", "-o", output_vcf
        ]

        logger.info(f"Computing GL for sample {name}, chromosome {chromosome}")
        mpileup_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        call_process = subprocess.run(call_command, stdin=mpileup_process.stdout, capture_output=True, text=True)
        mpileup_process.stdout.close()

        if call_process.returncode != 0:
            logger.error(f"Error in GL computation for sample {name}: {call_process.stderr}")
            raise RuntimeError(f"Error in GL computation for sample {name}: {call_process.stderr}")

        index_command = [TABIX, "-f", output_vcf]
        subprocess.run(index_command, check=True)


def merge_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    os.makedirs(glmergepath, exist_ok=True)

    gl_list_path = os.path.join(glpath, f"glimpse.{chromosome}_GL_list.txt")
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    with open(gl_list_path, "w") as gl_list:
        for file in os.listdir(glpath):
            if file.endswith(f".{chromosome}.vcf.gz"):
                gl_list.write(os.path.join(glpath, file) + "\n")

    command = [
        BCFTOOLS, "merge", "-@", f"{PARAMETERS['threads']}", "-m", "none", 
        "-r", chromosome, "-Oz", "-o", merged_vcf, "-l", gl_list_path
    ]

    logger.info(f"Merging GL files for chromosome {chromosome}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error merging GLs for chromosome {chromosome}: {process.stderr}")
        raise RuntimeError(f"Error merging GLs for chromosome {chromosome}: {process.stderr}")

    index_command = [TABIX, "-f", "-@", f"{PARAMETERS['threads']}", merged_vcf]
    subprocess.run(index_command, check=True)

    logger.info(f"Merged GL file created at {merged_vcf}")

def phase_genome(fq, chromosome):
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    os.makedirs(imputed_path, exist_ok=True)

    map_file = os.path.join(MAP_PATH, f"{chromosome}.b38.gmap.gz")
    reference_vcf = norm_vcf_path(chromosome)
    chunk_file = chunks_path(chromosome)
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    with open(chunk_file, "r") as chunks:
        for line in chunks:
            fields = line.strip().split()
            chunk_id = f"{int(fields[0]):02d}"
            input_region = fields[2]
            output_region = fields[3]
            output_vcf = os.path.join(imputed_path, f"glimpse.{chromosome}.{chunk_id}.imputed.vcf")

            command = [
                GLIMPSE_PHASE,
                "--input-gl", merged_vcf,
                "--reference", reference_vcf,
                "--map", map_file,
                "--input-region", input_region,
                "--output-region", output_region,
                "--output", output_vcf
            ]

            logger.info(f"Phasing chromosome {chromosome}, chunk {chunk_id}")
            process = subprocess.run(command, capture_output=True, text=True)
            if process.returncode != 0:
                logger.error(f"Error phasing chromosome {chromosome}, chunk {chunk_id}: {process.stderr}")
                raise RuntimeError(f"Error phasing chromosome {chromosome}, chunk {chunk_id}: {process.stderr}")

            bgzip_command = [BGZIP, output_vcf]
            tabix_command = [TABIX, "-f", "-@", f"{PARAMETERS['threads']}", f"{output_vcf}.gz"]

            subprocess.run(bgzip_command, check=True)
            subprocess.run(tabix_command, check=True)

def extract_chunk_id(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    merged_path = os.path.join(glimpse_outdir(fq), "imputed_file_merged")
    os.makedirs(merged_path, exist_ok=True)

    imputed_list = os.path.join(imputed_path, f"glimpse.{chromosome}_imputed_list.txt")

    # Hàm trích xuất chunk_id từ tên file
    def get_chunk_id(filename):
        match = re.search(rf"glimpse\.{chromosome}\.(\d+)\.imputed\.vcf\.gz", filename)
        return int(match.group(1)) if match else float('inf')

    # Lấy danh sách file phù hợp
    files = [
        file for file in os.listdir(imputed_path)
        if file.startswith(f"glimpse.{chromosome}.") and file.endswith(".imputed.vcf.gz")
    ]

    # Sắp xếp file theo chunk_id tăng dần
    sorted_files = sorted(files, key=get_chunk_id)

    # Ghi danh sách file đã sắp xếp vào imputed_list
    with open(imputed_list, "w") as imp_list:
        for file in sorted_files:
            imp_list.write(os.path.join(imputed_path, file) + "\n")


def ligate_genome(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    merged_path = os.path.join(glimpse_outdir(fq), "imputed_file_merged")
    os.makedirs(merged_path, exist_ok=True)

    imputed_list = os.path.join(imputed_path, f"glimpse.{chromosome}_imputed_list.txt")
    output_vcf = os.path.join(merged_path, f"glimpse.{chromosome}_imputed.vcf")

    command = [
        GLIMPSE_LIGATE, "--input", imputed_list, "--output", output_vcf
    ]

    logger.info(f"Ligating genome for chromosome {chromosome}")
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        logger.error(f"Error ligating chromosome {chromosome}: {process.stderr}")
        raise RuntimeError(f"Error ligating chromosome {chromosome}: {process.stderr}")

    bgzip_command = [BGZIP, output_vcf]
    tabix_command = [TABIX, "-f", "-@", f"{PARAMETERS['threads']}", f"{output_vcf}.gz"]

    subprocess.run(bgzip_command, check=True)
    subprocess.run(tabix_command, check=True)

def run_glimpse(fq):
    for chromosome in PARAMETERS["chrs"]:
        logger.info(f"Starting pipeline for chromosome {chromosome}...")

        # Step 1: Compute GLs
        compute_gls(fq, chromosome)

        # Step 2: Merge GLs
        merge_gls(fq, chromosome)

        # Step 3: Phase genome
        phase_genome(fq, chromosome)

        # Step 4: Ligate genome
        extract_chunk_id(fq, chromosome)
        ligate_genome(fq, chromosome)

        logger.info(f"Pipeline completed for chromosome {chromosome}.")
