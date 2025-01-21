import os
import subprocess
import shutil
from helper.config import TOOLS, PATHS, PARAMETERS
from helper.path_define import samid, tmp_outdir, batch1_final_outdir, bamlist_dir
from helper.logger import setup_logger

# Cấu hình từ JSON
REF = PATHS["ref"]
REF_INDEX_PREFIX = PATHS["ref_index_prefix"]
GATK_BUNDLE_DIR = PATHS["gatk_bundle_dir"]

logger = setup_logger(os.path.join(PATHS["logs"], "alignment_pipeline.log"))


def run_bwa_alignment(sample_id, fq, outdir, ref_index_prefix=REF, bwa=TOOLS["bwa"], samtools=TOOLS["samtools"]):
    """
    Runs a pipeline for alignment and BAM file processing using BWA and Samtools.
    """
    logger.info(f"Calculating {fq}. We'll save it in {outdir}")

    try:
        # Step 1: BWA alignment
        logger.info("\nRunning BWA alignment...")
        sai_file = os.path.join(outdir, f"{sample_id}.sai")
        bam_file = os.path.join(outdir, f"{sample_id}.bam")
        sorted_bam = os.path.join(outdir, f"{sample_id}.sorted.bam")
        rmdup_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.bam")
        finish_flag = os.path.join(outdir, "bwa_sort_rmdup.finish")

        bwa_aln_cmd = [
            bwa, "aln", "-e", "10", "-t", f"{PARAMETERS['threads']}", "-i", "5", "-q", "0",
            ref_index_prefix, fq
        ]
        with open(sai_file, "w") as sai_out:
            subprocess.run(bwa_aln_cmd, stdout=sai_out, check=True)

        bwa_samse_cmd = [
            bwa, "samse", "-r",
            f"@RG\\tID:default\\tPL:COMPLETE\\tSM:{sample_id}",
            ref_index_prefix, sai_file, fq
        ]

        samtools_view_cmd = [samtools, "view", "-h", "-Sb", "-@", f"{PARAMETERS['threads']}", "-"]

        # Open the output BAM file for writing
        with open(bam_file, "wb") as bam_out:
            # Start the BWA samse process
            bwa_process = subprocess.Popen(bwa_samse_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Start the Samtools view process, taking input from bwa_process
            subprocess.run(samtools_view_cmd, stdin=bwa_process.stdout, stdout=bam_out, check=True)
            
            # Close the stdout pipe to ensure proper cleanup
            bwa_process.stdout.close()
            bwa_process.wait()  # Ensure BWA process completes

        logger.info("** BWA done **")

        # Step 2: Sorting BAM
        logger.info("Sorting BAM...")
        subprocess.run([samtools, "sort", "-@", f"{PARAMETERS['threads']}", "-O", "bam", "-o", sorted_bam, bam_file], check=True)
        logger.info("** BAM sorted done **")

        # Step 3: Removing duplicates
        logger.info("Removing duplicates...")
        subprocess.run([samtools, "markdup", "-@", f"{PARAMETERS['threads']}", sorted_bam, rmdup_bam], check=True)
        logger.info("** rmdup done **")

        # Step 4: Indexing BAM
        logger.info("Indexing BAM...")
        subprocess.run([samtools, "index", "-@", f"{PARAMETERS['threads']}", rmdup_bam], check=True)
        logger.info("** index done **")

        # Step 5: Create finish flag
        with open(finish_flag, "w") as finish_file:
            finish_file.write("Pipeline completed successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"[WORKFLOW_ERROR_INFO] Command failed: {e.cmd}\nError: {e}")
        exit(1)

    # Check if the finish flag was created
    if not os.path.exists(finish_flag):
        logger.error("** [WORKFLOW_ERROR_INFO] bwa_sort_rmdup not done **")
        exit(1)

def run_bwa_realign(sample_id, outdir, ref=REF, gatk_bundle_dir=PATHS["gatk_bundle_dir"], gatk=TOOLS["gatk"], java=TOOLS["java"]):
    """
    Runs the RealignerTargetCreator step using GATK.
    """
    logger.info("\nStarting realign pipeline...")

    try:
        # Paths to input and output files
        bam_file = os.path.join(outdir, f"{sample_id}.sorted.rmdup.bam")
        intervals_file = os.path.join(outdir, f"{sample_id}.indel_target_intervals.list")
        realigned_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.bam")
        realigner_target_finish_flag = os.path.join(outdir, "RealignerTargetCreator.finish")
        indel_realigner_finish_flag = os.path.join(outdir, "IndelRealigner.finish")

        # Step 1: RealignerTargetCreator
        logger.info("Running RealignerTargetCreator...")
        realigner_target_cmd = [
            java, "-Xmx15g", "-jar", gatk,
            "-T", "RealignerTargetCreator",
            "-R", ref,
            "-I", bam_file,
            "-known", os.path.join(gatk_bundle_dir, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
            "-known", os.path.join(gatk_bundle_dir, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
            "-o", intervals_file
        ]
        subprocess.run(realigner_target_cmd, check=True)
        logger.info("** RealignerTargetCreator done **")

        # Create finish flag for RealignerTargetCreator
        with open(realigner_target_finish_flag, "w") as flag:
            flag.write("RealignerTargetCreator completed successfully.")

        # Verify the flag
        if not os.path.exists(realigner_target_finish_flag):
            raise FileNotFoundError("RealignerTargetCreator did not complete successfully.")

        # Step 2: IndelRealigner
        logger.info("Running IndelRealigner...")
        indel_realigner_cmd = [
            java, "-Xmx15g", "-jar", gatk,
            "-T", "IndelRealigner",
            "-R", ref,
            "-I", bam_file,
            "-known", os.path.join(gatk_bundle_dir, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
            "-known", os.path.join(gatk_bundle_dir, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
            "--targetIntervals", intervals_file,
            "-o", realigned_bam
        ]
        subprocess.run(indel_realigner_cmd, check=True)
        logger.info("** IndelRealigner done **")

        # Create finish flag for IndelRealigner
        with open(indel_realigner_finish_flag, "w") as flag:
            flag.write("IndelRealigner completed successfully.")

        # Verify the flag
        if not os.path.exists(indel_realigner_finish_flag):
            raise FileNotFoundError("IndelRealigner did not complete successfully.")

    except subprocess.CalledProcessError as e:
        logger.info(f"[WORKFLOW_ERROR_INFO] Command failed: {e.cmd}\nError: {e}")
        exit(1)
    except FileNotFoundError as e:
        logger.info(f"[WORKFLOW_ERROR_INFO] {e}")
        exit(1)

    logger.info("Realign pipeline completed successfully.")

def run_bqsr(sample_id, outdir, ref=REF, gatk_bundle_dir=PATHS["gatk_bundle_dir"], gatk=TOOLS["gatk"], samtools=TOOLS["samtools"], java=TOOLS["java"]):
    """
    Runs the Base Quality Score Recalibration (BQSR) pipeline using GATK and Samtools.

    Parameters:
        sample_id (str): Sample identifier.
        ref (str): Path to the reference genome.
        gatk_bundle_dir (str): Path to the GATK bundle directory.
        outdir (str): Directory to save the output files.
        gatk (str): Path to the GATK executable (default: "gatk").
        samtools (str): Path to the Samtools executable (default: "samtools").
        java (str): Path to the Java executable (default: "java").
    """
    try:
        # Paths to files
        realigned_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.bam")
        recal_table = os.path.join(outdir, f"{sample_id}.recal_data.table")
        bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
        index_flag = os.path.join(outdir, "index_sorted_rmdup_realign_bam.finish")
        recal_flag = os.path.join(outdir, "baseRecalibrator.finish")
        print_reads_flag = os.path.join(outdir, "PrintReads.finish")
        bam_index_flag = os.path.join(outdir, "bam_index.finish")

        # Step 1: Index the realigned BAM
        logger.info("Indexing realigned BAM...")
        subprocess.run([samtools, "index", "-@", f"{PARAMETERS['threads']}", realigned_bam], check=True)
        logger.info("** Index done **")
        with open(index_flag, "w") as flag:
            flag.write("Indexing completed successfully.")

        if not os.path.exists(index_flag):
            raise FileNotFoundError("Indexing of realigned BAM did not complete successfully.")

        # Step 2: BaseRecalibrator
        logger.info("Running BaseRecalibrator...")
        base_recal_cmd = [
            java, "-jar", gatk,
            "-T", "BaseRecalibrator",
            "-nct", "8",
            "-R", ref,
            "-I", realigned_bam,
            "--knownSites", os.path.join(gatk_bundle_dir, "Homo_sapiens_assembly38.dbsnp138.vcf.gz"),
            "--knownSites", os.path.join(gatk_bundle_dir, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
            "--knownSites", os.path.join(gatk_bundle_dir, "Homo_sapiens_assembly38.known_indels.vcf.gz"),
            "-o", recal_table
        ]
        subprocess.run(base_recal_cmd, check=True)
        logger.info("** BaseRecalibrator done **")
        with open(recal_flag, "w") as flag:
            flag.write("BaseRecalibrator completed successfully.")

        if not os.path.exists(recal_flag):
            raise FileNotFoundError("BaseRecalibrator did not complete successfully.")

        # Step 3: PrintReads
        logger.info("Running PrintReads...")
        print_reads_cmd = [
            java, "-jar", gatk,
            "-T", "PrintReads",
            "-nct", "8",
            "-R", ref,
            "--BQSR", recal_table,
            "-I", realigned_bam,
            "-o", bqsr_bam
        ]
        subprocess.run(print_reads_cmd, check=True)
        logger.info("** PrintReads done **")
        with open(print_reads_flag, "w") as flag:
            flag.write("PrintReads completed successfully.")

        if not os.path.exists(print_reads_flag):
            raise FileNotFoundError("PrintReads did not complete successfully.")

        # Step 4: Index the BQSR BAM
        logger.info("Indexing BQSR BAM...")
        subprocess.run([samtools, "index", "-@", f"{PARAMETERS['threads']}", bqsr_bam], check=True)
        logger.info("** BAM index done **")
        with open(bam_index_flag, "w") as flag:
            flag.write("BAM indexing completed successfully.")

        if not os.path.exists(bam_index_flag):
            raise FileNotFoundError("BAM indexing of BQSR BAM did not complete successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"[WORKFLOW_ERROR_INFO] Command failed: {e.cmd}\nError: {e}")
        exit(1)
    except FileNotFoundError as e:
        logger.error(f"[WORKFLOW_ERROR_INFO] {e}")
        exit(1)

    logger.info("BQSR pipeline completed successfully.")

def run_bam_stats(sample_id, outdir, samtools=TOOLS["samtools"]):
    """
    Runs Samtools stats on the BQSR BAM file and generates statistics.
    """
    try:
        # Paths to input and output files
        bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
        bam_stats_file = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bamstats")
        bam_stats_flag = os.path.join(outdir, "bamstats.finish")

        # Step: Run Samtools stats
        logger.info("Running Samtools stats...")
        with open(bam_stats_file, "w") as stats_out:
            subprocess.run([samtools, "stats", "-@", f"{PARAMETERS['threads']}", bqsr_bam], stdout=stats_out, check=True)
        logger.info("** bamstats done **")

        # Create finish flag
        with open(bam_stats_flag, "w") as flag:
            flag.write("bamstats completed successfully.")

        # Verify the flag
        if not os.path.exists(bam_stats_flag):
            raise FileNotFoundError("bamstats did not complete successfully.")

    except subprocess.CalledProcessError as e:
        logger.info(f"[WORKFLOW_ERROR_INFO] Command failed: {e.cmd}\nError: {e}")
        exit(1)
    except FileNotFoundError as e:
        logger.info(f"[WORKFLOW_ERROR_INFO] {e}")
        exit(1)

    logger.info("BAM stats pipeline completed successfully.")


def run_bedtools(fq, sample_id, outdir, final_outdir, bedtools=TOOLS["bedtools"], bgzip=TOOLS["bgzip"], tabix=TOOLS["tabix"]):
    """
    Runs Bedtools genome coverage analysis and compresses the output using bgzip.
    """
    try:
        # Paths to input and output files
        bqsr_bam = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.bam")
        cvg_bed_gz = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz")
        finish_flag = os.path.join(outdir, "sorted_rmdup_realign_BQSR_cvg_bed.finish")
        bam_list_file = bamlist_dir(fq)

        # Step 1: Bedtools genome coverage
        logger.info("Running Bedtools genome coverage...")
        bedtools_cmd = [
            bedtools, "genomecov", "-ibam", bqsr_bam, "-bga", "-split"
        ]
        with open(cvg_bed_gz, "wb") as cvg_out:
            subprocess.run(bedtools_cmd, stdout=subprocess.PIPE, check=True)
            subprocess.run([bgzip], stdin=subprocess.PIPE, stdout=cvg_out, check=True)
        logger.info("** sorted.rmdup.realign.BQSR.cvg.bed.gz done **")

        # Step 2: Index the compressed BED file
        logger.info("Indexing the compressed BED file with Tabix...")
        subprocess.run([tabix, "-p", "bed", "-@", f"{PARAMETERS['threads']}", cvg_bed_gz], check=True)

        # Create finish flag
        with open(finish_flag, "w") as flag:
            flag.write("Bedtools pipeline completed successfully.")

        if not os.path.exists(finish_flag):
            raise FileNotFoundError("Bedtools pipeline did not complete successfully.")

        # Step 3: Move files to the final output directory, generate bam.list 
        logger.info("Moving final files to the output directory...")
        for file_suffix in [".bam", ".bam.bai", ".cvg.bed.gz", ".cvg.bed.gz.tbi"]:
            src_file = os.path.join(outdir, f"{sample_id}.sorted.rmdup.realign.BQSR{file_suffix}")
            dst_file = os.path.join(final_outdir, os.path.basename(src_file))
            if os.path.exists(src_file):
                os.rename(src_file, dst_file)
                if file_suffix == ".bam":
                    with open(bam_list_file, "a") as bam_list:
                        bam_list.write(f"{dst_file}\n")

        # Step 4: Remove the temporary output directory
        logger.info("Removing temporary output directory...")
        shutil.rmtree(outdir)
        logger.info(f"Temporary directory {outdir} deleted.")

    except subprocess.CalledProcessError as e:
        logger.error(f"[WORKFLOW_ERROR_INFO] Command failed: {e.cmd}\nError: {e}")
        exit(1)
    except FileNotFoundError as e:
        logger.error(f"[WORKFLOW_ERROR_INFO] {e}")
        exit(1)

    logger.info("Bedtools pipeline completed successfully.")


def run_alignment_pipeline(fq):
    """
    Thực hiện pipeline alignment cho một mẫu FASTQ.
    """

    if os.path.exists(batch1_final_outdir(fq)):
        logger.info(f"Đã có thư mục kết quả alignment cho mẫu {samid(fq)}")
        return

    # Tạo các thư mục nếu chưa tồn tại
    os.makedirs(tmp_outdir(fq), exist_ok=True)
    os.makedirs(batch1_final_outdir(fq), exist_ok=True)

    logger.info(f"=== Bắt đầu pipeline alignment cho mẫu {samid(fq)} ===")
    logger.info(f"FASTQ: {fq}")
    logger.info(f"Thư mục tạm: {tmp_outdir(fq)}")
    logger.info(f"Thư mục kết quả cuối: {batch1_final_outdir(fq)}")

    # Step 1: Chạy BWA để căn chỉnh và loại bỏ bản sao (duplicates)
    run_bwa_alignment(samid(fq), fq, tmp_outdir(fq))
    logger.info(f"Hoàn thành BWA alignment. Sample: {samid(fq)}")

    # Step 2: Thực hiện realignment
    run_bwa_realign(samid(fq), tmp_outdir(fq))
    logger.info(f"Hoàn thành tmp_outdir(fq). Sample: {samid(fq)}")

    # Step 3: Recalibrate Base Quality Scores (BQSR)
    run_bqsr(samid(fq), tmp_outdir(fq))
    logger.info(f"Hoàn thành BQSR. Sample: {samid(fq)}")

    # Step 4: Tạo thống kê BAM và coverage
    run_bam_stats(samid(fq), tmp_outdir(fq))
    logger.info(f"Hoàn thành thống kê và coverage cho BAM.")

    # Step 5: Di chuyển file kết quả cuối cùng vào batch1_final_files
    run_bedtools(fq, samid(fq), tmp_outdir(fq), batch1_final_outdir(fq))
    logger.info(f"Kết quả đã được lưu tại {batch1_final_outdir(fq)}")

    logger.info(f"=== Hoàn thành pipeline alignment cho mẫu {samid(fq)} ===")
