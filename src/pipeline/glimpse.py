import os
from helper.slurm import create_slurm_script, submit_to_slurm, wait_for_job_completion
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import bamlist_dir, glimpse_outdir
from helper.path_define import filtered_vcf_path, filtered_tsv_path, chunks_path

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

        job_name = f"GL_{name}"
        script_name = f"glimpse.s2.gl.{name}.{chromosome}.slurm"
        output_vcf = os.path.join(glpath, f"{name}.{chromosome}.vcf.gz")

        commands = f"""
VCF={filtered_vcf_path(chromosome)}
TSV={filtered_tsv_path(chromosome)}

{BCFTOOLS} mpileup -f {REF} -I -E -a 'FORMAT/DP' -T $VCF -r {chromosome} {line} -Ou | \
    {BCFTOOLS} call -Aim -C alleles -T $TSV -Oz -o {output_vcf}
{BCFTOOLS} index -f {output_vcf}
        """

        create_slurm_script(
            script_name=script_name,
            job_name=job_name,
            commands=commands,
            output_dir=glpath,
            cpus_per_task=2,
            mem=20480,
            time="7-00:00:00",
        )
        job_id = submit_to_slurm(script_name)
        print(f"Submitted GL computation job for sample {name}, chromosome {chromosome} with Job ID: {job_id}")
        return job_id

def merge_gls(fq, chromosome):
    glpath = os.path.join(glimpse_outdir(fq), "GL_file")
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    os.makedirs(glmergepath, exist_ok=True)

    job_name = f"GL_merge_{chromosome}"
    script_name = f"glimpse.s3.GL_merge.{chromosome}.slurm"
    gl_list_path = os.path.join(glpath, f"glimpse.{chromosome}_GL_list.txt")
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    commands = f"""
ls {glpath}/*.{chromosome}.vcf.gz > {gl_list_path}
{BCFTOOLS} merge -m none -r {chromosome} -Oz -o {merged_vcf} -l {gl_list_path}
{BCFTOOLS} index -f {merged_vcf}
    """

    create_slurm_script(
        script_name=script_name,
        job_name=job_name,
        commands=commands,
        output_dir=glmergepath,
        cpus_per_task=4,
        mem=20480,
        time="7-00:00:00",
    )
    job_id = submit_to_slurm(script_name)
    print(f"Submitted GL merge job for chromosome {chromosome} with Job ID: {job_id}")
    return job_id

def phase_genome(fq, chromosome):
    glmergepath = os.path.join(glimpse_outdir(fq), "GL_file_merged")
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    os.makedirs(imputed_path, exist_ok=True)

    map_file = os.path.join(MAP_PATH, f"{chromosome}.b38.gmap.gz")
    reference_vcf = filtered_vcf_path(chromosome)
    chunk_file = chunks_path(chromosome)
    merged_vcf = os.path.join(glmergepath, f"glimpse.{chromosome}.vcf.gz")

    with open(chunk_file, "r") as chunks:
        for line in chunks:
            fields = line.strip().split()
            chunk_id = f"{int(fields[0]):02d}"
            input_region = fields[2]
            output_region = fields[3]
            output_vcf = os.path.join(imputed_path, f"glimpse.{chromosome}.{chunk_id}.imputed.vcf")

            job_name = f"phase_chr{chromosome}_{chunk_id}"
            script_name = f"glimpse.s4.phase.{chromosome}.{chunk_id}.slurm"

            commands = f"""
{GLIMPSE_PHASE} --input-gl {merged_vcf} --reference {reference_vcf} --map {map_file} \
    --input-region {input_region} --output-region {output_region} --output {output_vcf}
{BGZIP} {output_vcf}
{BCFTOOLS} index -f {output_vcf}.gz
            """

            create_slurm_script(
                script_name=script_name,
                job_name=job_name,
                commands=commands,
                output_dir=imputed_path,
                cpus_per_task=1,
                mem=20480,
                time="7-00:00:00",
            )
            job_id = submit_to_slurm(script_name)
            print(f"Submitted phasing job for chromosome {chromosome}, chunk {chunk_id} with Job ID: {job_id}")
            return job_id

def ligate_genome(fq, chromosome):
    imputed_path = os.path.join(glimpse_outdir(fq), "imputed_file")
    merged_path = os.path.join(glimpse_outdir(fq), "imputed_file_merged")
    os.makedirs(merged_path, exist_ok=True)

    imputed_list = os.path.join(imputed_path, f"glimpse.{chromosome}_imputed_list.txt")
    output_vcf = os.path.join(merged_path, f"glimpse.{chromosome}_imputed.vcf")

    job_name = f"ligate_chr{chromosome}"
    script_name = f"glimpse.s5.ligate.{chromosome}.slurm"

    commands = f"""
ls {imputed_path}/glimpse.{chromosome}.*.imputed.vcf.gz > {imputed_list}
{GLIMPSE_LIGATE} --input {imputed_list} --output {output_vcf}
{BGZIP} {output_vcf}
{BCFTOOLS} index -f {output_vcf}.gz
    """

    create_slurm_script(
        script_name=script_name,
        job_name=job_name,
        commands=commands,
        output_dir=merged_path,
        cpus_per_task=8,
        mem=20480,
        time="7-00:00:00",
    )
    job_id = submit_to_slurm(script_name)
    print(f"Submitted ligation job for chromosome {chromosome} with Job ID: {job_id}")
    return job_id

def run_pipeline(fq):
    for chromosome in PARAMETERS["chrs"]:
        print(f"Starting pipeline for chromosome {chromosome}...")

        # Step 1: Compute GLs
        gl_job_id = compute_gls(fq, chromosome)
        wait_for_job_completion(gl_job_id)

        # Step 2: Merge GLs
        merge_job_id = merge_gls(fq, chromosome)
        wait_for_job_completion(merge_job_id)

        # Step 3: Phase genome
        phase_job_id = phase_genome(fq, chromosome)
        wait_for_job_completion(phase_job_id)

        # Step 4: Ligate genome
        ligate_job_id = ligate_genome(fq, chromosome)
        wait_for_job_completion(ligate_job_id)

        print(f"Pipeline completed for chromosome {chromosome}.")
