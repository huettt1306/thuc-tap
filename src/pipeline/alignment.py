import os, subprocess
from helper.config import TOOLS, PATHS
from helper.path_define import samid, tmp_outdir, batch1_final_outdir
from helper.logger import setup_logger

# Cấu hình từ JSON
REF = PATHS["ref"]
REF_INDEX_PREFIX = PATHS["ref_index_prefix"]
GATK_BUNDLE_DIR = PATHS["gatk_bundle_dir"]

logger = setup_logger(log_file="logs/alignment.log")

def run_bwa_alignment(sample_id, fq, outdir):
    """
    Chạy BWA alignment pipeline từ FASTQ sang BAM.
    """
    logger.info(f"Running BWA alignment for {sample_id}...")

    sai_file = f"{outdir}/{sample_id}.sai"
    bam_file = f"{outdir}/{sample_id}.bam"
    sorted_bam = f"{outdir}/{sample_id}.sorted.bam"
    rmdup_bam = f"{outdir}/{sample_id}.sorted.rmdup.bam"

    try:
        # BWA alignment
        bwa_aln_command = [
            TOOLS["bwa"], "aln", "-e", "10", "-t", "8", "-i", "5", "-q", "0",
            PATHS["ref_index_prefix"], fq
        ]
        with open(sai_file, "w") as sai_output:
            subprocess.run(bwa_aln_command, stdout=sai_output, check=True, text=True)
        logger.info(f"BWA alignment completed for {sample_id}, output: {sai_file}")

        # BWA samse and BAM conversion
        bwa_samse_command = [
            TOOLS["bwa"], "samse", "-r", f"@RG\tID:{sample_id}\tPL:COMPLETE\tSM:{sample_id}",
            PATHS["ref_index_prefix"], sai_file, fq
        ]
        samtools_view_command = [
            TOOLS["samtools"], "view", "-h", "-Sb", "-o", bam_file, "-"
        ]
        bwa_process = subprocess.Popen(bwa_samse_command, stdout=subprocess.PIPE, text=True)
        subprocess.run(samtools_view_command, stdin=bwa_process.stdout, check=True, text=True)
        logger.info(f"BAM file created: {bam_file}")

        # Sort BAM
        samtools_sort_command = [
            TOOLS["samtools"], "sort", "-@", "8", "-O", "bam", "-o", sorted_bam, bam_file
        ]
        subprocess.run(samtools_sort_command, check=True, text=True)
        logger.info(f"Sorted BAM file created: {sorted_bam}")

        # Remove duplicates
        samtools_rmdup_command = [
            TOOLS["samtools"], "rmdup", sorted_bam, rmdup_bam
        ]
        subprocess.run(samtools_rmdup_command, check=True, text=True)
        logger.info(f"Removed duplicates, output: {rmdup_bam}")

        # Index BAM
        samtools_index_command = [
            TOOLS["samtools"], "index", rmdup_bam
        ]
        subprocess.run(samtools_index_command, check=True, text=True)
        logger.info(f"Index created for BAM file: {rmdup_bam}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error during BWA alignment pipeline for {sample_id}: {e}")
        raise


def run_realign_pipeline(sample_id, outdir):
    """
    Chạy GATK realignment pipeline từ BAM đầu vào.
    """
    logger.info(f"Running GATK realignment for {sample_id}...")

    bam_file = f"{outdir}/{sample_id}.sorted.rmdup.bam"
    intervals_file = f"{outdir}/{sample_id}.indel_target_intervals.list"
    realigned_bam = f"{outdir}/{sample_id}.sorted.rmdup.realign.bam"

    try:
        # RealignerTargetCreator
        target_creator_command = [
            TOOLS["java"], "-Xmx15g", "-jar", TOOLS["gatk"],
            "-T", "RealignerTargetCreator",
            "-R", PATHS["ref"],
            "-I", bam_file,
            "-known", f"{PATHS['gatk_bundle']}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "-known", f"{PATHS['gatk_bundle']}/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "-o", intervals_file
        ]
        subprocess.run(target_creator_command, check=True, text=True)
        logger.info(f"Indel target intervals created: {intervals_file}")

        # IndelRealigner
        realigner_command = [
            TOOLS["java"], "-Xmx15g", "-jar", TOOLS["gatk"],
            "-T", "IndelRealigner",
            "-R", PATHS["ref"],
            "-I", bam_file,
            "-known", f"{PATHS['gatk_bundle']}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "-known", f"{PATHS['gatk_bundle']}/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "--targetIntervals", intervals_file,
            "-o", realigned_bam
        ]
        subprocess.run(realigner_command, check=True, text=True)
        logger.info(f"Realigned BAM file created: {realigned_bam}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error during GATK realignment pipeline for {sample_id}: {e}")
        raise


def run_bqsr_pipeline(sample_id, outdir):
    """
    Chạy GATK Base Quality Score Recalibration (BQSR) pipeline.
    """
    logger.info(f"Running BQSR pipeline for {sample_id}...")

    realigned_bam = f"{outdir}/{sample_id}.sorted.rmdup.realign.bam"
    recal_table = f"{outdir}/{sample_id}.recal_data.table"
    bqsr_bam = f"{outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bam"

    try:
        # BaseRecalibrator
        base_recalibrator_command = [
            TOOLS["java"], "-jar", TOOLS["gatk"],
            "-T", "BaseRecalibrator",
            "-nct", "8",
            "-R", PATHS["ref"],
            "-I", realigned_bam,
            "--knownSites", f"{PATHS['gatk_bundle']}/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
            "--knownSites", f"{PATHS['gatk_bundle']}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "--knownSites", f"{PATHS['gatk_bundle']}/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "-o", recal_table
        ]
        subprocess.run(base_recalibrator_command, check=True, text=True)
        logger.info(f"Recalibration table created: {recal_table}")

        # PrintReads
        print_reads_command = [
            TOOLS["java"], "-jar", TOOLS["gatk"],
            "-T", "PrintReads",
            "-nct", "8",
            "-R", PATHS["ref"],
            "--BQSR", recal_table,
            "-I", realigned_bam,
            "-o", bqsr_bam
        ]
        subprocess.run(print_reads_command, check=True, text=True)
        logger.info(f"BQSR BAM file created: {bqsr_bam}")

        # Index BAM
        samtools_index_command = [
            TOOLS["samtools"], "index", bqsr_bam
        ]
        subprocess.run(samtools_index_command, check=True, text=True)
        logger.info(f"Index created for BQSR BAM file: {bqsr_bam}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error during BQSR pipeline for {sample_id}: {e}")
        raise


def run_bam_stats_pipeline(sample_id, outdir):
    """
    Chạy BAM stats pipeline để tạo thống kê BAM.
    """
    logger.info(f"Running BAM stats for {sample_id}...")

    bam_file = f"{outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bam"
    stats_file = f"{outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bamstats"

    try:
        # BAM stats
        samtools_stats_command = [
            TOOLS["samtools"], "stats", bam_file
        ]
        with open(stats_file, "w") as stats_output:
            subprocess.run(samtools_stats_command, stdout=stats_output, check=True, text=True)
        logger.info(f"BAM stats created: {stats_file}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error during BAM stats pipeline for {sample_id}: {e}")
        raise

def move_final_files_pipeline(sample_id, tmp_dir, final_dir):
    """
    Di chuyển các file kết quả cuối cùng từ thư mục tạm sang thư mục đích.
    """
    logger.info(f"Moving final files for {sample_id}...")

    try:
        for ext in ["bam", "bam.bai", "bamstats"]:
            src_file = f"{tmp_dir}/{sample_id}.sorted.rmdup.realign.BQSR.{ext}"
            dest_file = f"{final_dir}/{sample_id}.sorted.rmdup.realign.BQSR.{ext}"
            if os.path.exists(src_file):
                subprocess.run(["mv", src_file, dest_file], check=True)
                logger.info(f"Moved {src_file} to {dest_file}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error moving final files for {sample_id}: {e}")
        raise


def run_alignment_pipeline(fq):
    """
    Thực hiện pipeline alignment cho một mẫu FASTQ.
    """
    sample_id = samid(fq)
    tmp_dir = tmp_outdir(fq)
    final_dir = batch1_final_outdir(fq)

    # Tạo các thư mục nếu chưa tồn tại
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)

    logger.info(f"=== Bắt đầu pipeline alignment cho mẫu {sample_id} ===")
    logger.info(f"FASTQ: {fq}")
    logger.info(f"Thư mục tạm: {tmp_dir}")
    logger.info(f"Thư mục kết quả cuối: {final_dir}")

    try:
        # Step 1: Submit BWA alignment
        run_bwa_alignment(sample_id, fq, tmp_dir)

        # Step 2: Submit GATK realignment
        run_realign_pipeline(sample_id, tmp_dir)

        # Step 3: Submit BQSR
        run_bqsr_pipeline(sample_id, tmp_dir)

        # Step 4: Submit BAM stats
        run_bam_stats_pipeline(sample_id, tmp_dir)

        # Step 5: Submit file move job
        move_final_files_pipeline(sample_id, tmp_dir, final_dir)

        logger.info(f"=== Hoàn thành pipeline alignment cho mẫu {sample_id} ===")
    except Exception as e:
        logger.error(f"Error during alignment pipeline for {sample_id}: {e}")
        raise
