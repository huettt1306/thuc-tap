import os
from helper.config import TOOLS, PATHS
from helper.path_define import samid, tmp_outdir, batch1_final_outdir
from helper.slurm import create_slurm_script, submit_to_slurm, wait_for_job_completion

# Cấu hình từ JSON
REF = PATHS["ref"]
REF_INDEX_PREFIX = PATHS["ref_index_prefix"]
GATK_BUNDLE_DIR = PATHS["gatk_bundle_dir"]

def run_bwa_alignment(sample_id, fq, outdir):
    """
    Submit BWA alignment job to SLURM.
    """
    print(f"Preparing SLURM script for BWA alignment of {fq}...")
    commands = f"""
sai_file={outdir}/{sample_id}.sai
bam_file={outdir}/{sample_id}.bam
sorted_bam={outdir}/{sample_id}.sorted.bam
rmdup_bam={outdir}/{sample_id}.sorted.rmdup.bam

# BWA alignment
{TOOLS['bwa']} aln -e 10 -t 8 -i 5 -q 0 {REF_INDEX_PREFIX} {fq} > $sai_file

# BWA samse and BAM conversion
{TOOLS['bwa']} samse -r "@RG\\tID:{sample_id}\\tPL:COMPLETE\\tSM:{sample_id}" {REF_INDEX_PREFIX} $sai_file {fq} | \
    {TOOLS['samtools']} view -h -Sb - > $bam_file

# Sort BAM
{TOOLS['samtools']} sort -@ 8 -O bam -o $sorted_bam $bam_file

# Remove duplicates
{TOOLS['samtools']} rmdup $sorted_bam $rmdup_bam

# Index BAM
{TOOLS['samtools']} index $rmdup_bam
"""
    script_name = f"{outdir}/{sample_id}_bwa_alignment.slurm"
    job_name = f"{sample_id}_bwa_alignment"

    # Create and submit SLURM script
    create_slurm_script(script_name, job_name, commands, outdir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted BWA alignment job for {sample_id} with Job ID: {job_id}")
    return job_id

def run_realign_pipeline(sample_id, outdir):
    """
    Submit GATK Realignment job to SLURM.
    """
    print(f"Preparing SLURM script for GATK realignment of {sample_id}...")
    commands = f"""
bam_file={outdir}/{sample_id}.sorted.rmdup.bam
intervals_file={outdir}/{sample_id}.indel_target_intervals.list
realigned_bam={outdir}/{sample_id}.sorted.rmdup.realign.bam

# RealignerTargetCreator
{TOOLS['java']} -Xmx15g -jar {TOOLS['gatk']} \
    -T RealignerTargetCreator \
    -R {REF} \
    -I $bam_file \
    -known {GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known {GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $intervals_file

# IndelRealigner
{TOOLS['java']} -Xmx15g -jar {TOOLS['gatk']} \
    -T IndelRealigner \
    -R {REF} \
    -I $bam_file \
    -known {GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known {GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --targetIntervals $intervals_file \
    -o $realigned_bam
"""
    script_name = f"{outdir}/{sample_id}_realign.slurm"
    job_name = f"{sample_id}_realign"

    # Create and submit SLURM script
    create_slurm_script(script_name, job_name, commands, outdir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted GATK realignment job for {sample_id} with Job ID: {job_id}")
    return job_id

def run_bqsr_pipeline(sample_id, outdir):
    """
    Submit Base Quality Score Recalibration (BQSR) job to SLURM.
    """
    print(f"Preparing SLURM script for BQSR of {sample_id}...")
    commands = f"""
realigned_bam={outdir}/{sample_id}.sorted.rmdup.realign.bam
recal_table={outdir}/{sample_id}.recal_data.table
bqsr_bam={outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bam

# BaseRecalibrator
{TOOLS['java']} -jar {TOOLS['gatk']} \
    -T BaseRecalibrator \
    -nct 8 \
    -R {REF} \
    -I $realigned_bam \
    --knownSites {GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --knownSites {GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites {GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $recal_table

# PrintReads
{TOOLS['java']} -jar {TOOLS['gatk']} \
    -T PrintReads \
    -nct 8 \
    -R {REF} \
    --BQSR $recal_table \
    -I $realigned_bam \
    -o $bqsr_bam

# Index BAM
{TOOLS['samtools']} index $bqsr_bam
"""
    script_name = f"{outdir}/{sample_id}_bqsr.slurm"
    job_name = f"{sample_id}_bqsr"

    # Create and submit SLURM script
    create_slurm_script(script_name, job_name, commands, outdir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted BQSR job for {sample_id} with Job ID: {job_id}")
    return job_id

def run_bam_stats_pipeline(sample_id, outdir):
    """
    Submit BAM stats job to SLURM.
    """
    print(f"Preparing SLURM script for BAM stats of {sample_id}...")
    commands = f"""
bam_file={outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bam
stats_file={outdir}/{sample_id}.sorted.rmdup.realign.BQSR.bamstats

# BAM stats
{TOOLS['samtools']} stats $bam_file > $stats_file
"""
    script_name = f"{outdir}/{sample_id}_bam_stats.slurm"
    job_name = f"{sample_id}_bam_stats"

    # Create and submit SLURM script
    create_slurm_script(script_name, job_name, commands, outdir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted BAM stats job for {sample_id} with Job ID: {job_id}")
    return job_id

def move_final_files_pipeline(sample_id, tmp_dir, final_dir):
    """
    Submit file move job to SLURM.
    """
    print(f"Preparing SLURM script to move final files for {sample_id}...")
    commands = f"""
# Move final files
for ext in bam bam.bai bamstats; do
    if [ -f {tmp_dir}/{sample_id}.sorted.rmdup.realign.BQSR.$ext ]; then
        mv {tmp_dir}/{sample_id}.sorted.rmdup.realign.BQSR.$ext {final_dir}/
    fi
done
"""
    script_name = f"{tmp_dir}/{sample_id}_move_files.slurm"
    job_name = f"{sample_id}_move_files"

    # Create and submit SLURM script
    create_slurm_script(script_name, job_name, commands, tmp_dir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted file move job for {sample_id} with Job ID: {job_id}")
    return job_id

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

    print(f"=== Bắt đầu pipeline alignment cho mẫu {sample_id} ===")
    print(f"FASTQ: {fq}")
    print(f"Thư mục tạm: {tmp_dir}")
    print(f"Thư mục kết quả cuối: {final_dir}")

    # Step 1: Submit BWA alignment
    bwa_job_id = run_bwa_alignment(sample_id, fq, tmp_dir)
    wait_for_job_completion(bwa_job_id)

    # Step 2: Submit GATK realignment
    realign_job_id = run_realign_pipeline(sample_id, tmp_dir)
    wait_for_job_completion(realign_job_id)

    # Step 3: Submit BQSR
    bqsr_job_id = run_bqsr_pipeline(sample_id, tmp_dir)
    wait_for_job_completion(bqsr_job_id)

    # Step 4: Submit BAM stats
    bam_stats_job_id = run_bam_stats_pipeline(sample_id, tmp_dir)
    wait_for_job_completion(bam_stats_job_id)

    # Step 5: Submit file move job
    move_files_job_id = move_final_files_pipeline(sample_id, tmp_dir, final_dir)
    wait_for_job_completion(move_files_job_id)

    print(f"=== Hoàn thành pipeline alignment cho mẫu {sample_id} ===")
