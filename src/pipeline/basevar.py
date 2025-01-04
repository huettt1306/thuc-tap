import os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.path_define import basevar_outdir, bamlist_dir, vcf_list_path, basevar_vcf
from helper.slurm import create_slurm_script, submit_to_slurm, wait_for_job_completion

# Load global configurations
REF = PATHS["ref"]
REF_FAI = PATHS["ref_fai"]
BCFTOOLS = TOOLS["bcftools"]
TABIX = TOOLS["tabix"]
DELTA = PARAMETERS["basevar"]["delta"]
THREADS = PARAMETERS["basevar"]["threads"]


def load_reference_fai(in_fai, chroms=None):
    ref = []
    with open(in_fai) as fh:
        for r in fh:
            col = r.strip().split()
            if chroms is not None and len(chroms):
                if col[0] in chroms:
                    ref.append([col[0], 1, int(col[1])])
            else:
                ref.append([col[0], 1, int(col[1])])
    return ref


def run_basevar_step(fq, chromosome):
    ref_fai = load_reference_fai(REF_FAI, [chromosome])
    bamlist_path = bamlist_dir(fq)
    outdir = basevar_outdir(fq)

    jobs = []
    for chr_id, reg_start, reg_end in ref_fai:
        for i in range(reg_start - 1, reg_end, DELTA):
            start = i + 1
            end = i + DELTA if i + DELTA <= reg_end else reg_end
            region = f"{chr_id}:{start}-{end}"
            outfile_prefix = f"{chr_id}_{start}_{end}"

            commands = f"""
{TOOLS['basevar']} basetype \
    -t {THREADS} \
    -R {REF} \
    -L {bamlist_path} \
    -r {region} \
    --min-af 0.001 \
    --output-vcf {outdir}/{outfile_prefix}.vcf.gz \
    --output-cvg {outdir}/{outfile_prefix}.cvg.tsv.gz \
    --smart-rerun > {outdir}/{outfile_prefix}.log && \
    echo "** {outfile_prefix} done **"
"""

            script_name = os.path.join(outdir, f"basevar_{outfile_prefix}.slurm")
            job_name = f"basevar_{outfile_prefix}"

            create_slurm_script(script_name, job_name, commands, outdir)
            job_id = submit_to_slurm(script_name)
            print(f"Submitted BaseVar job for {region} with Job ID: {job_id}")
            jobs.append(job_id)

    return jobs


def create_vcf_list(fq, chromosome):
    vcf_list = vcf_list_path(fq, chromosome)
    print(f"=== Creating VCF list for chromosome {chromosome} ===")
    with open(vcf_list, "w") as f:
        for root, _, files in os.walk(basevar_outdir(fq)):
            for file in files:
                if file.endswith(".vcf.gz") and f"{chromosome}_" in file:
                    f.write(os.path.join(root, file) + "\n")
    print(f"VCF list saved at {vcf_list}")
    return vcf_list


def merge_vcf_files_slurm(fq, chromosome):
    outdir = basevar_outdir(fq)
    vcf_list_path = create_vcf_list(fq, chromosome)
    merged_vcf = basevar_vcf(fq, chromosome)

    print(f"Preparing SLURM script for merging VCF files for chromosome {chromosome}...")
    commands = f"""
{BCFTOOLS} concat \
    --threads {THREADS} \
    -a --rm-dups all \
    -O z \
    -o {merged_vcf} \
    $(cat {vcf_list_path})
"""

    script_name = os.path.join(outdir, f"merge_vcf_{chromosome}.slurm")
    job_name = f"merge_vcf_{chromosome}"

    create_slurm_script(script_name, job_name, commands, outdir)
    job_id = submit_to_slurm(script_name)
    print(f"Submitted VCF merge job for chromosome {chromosome} with Job ID: {job_id}")
    return job_id, merged_vcf


def index_vcf_file_slurm(vcf_path):
    print(f"Preparing SLURM script for indexing VCF file {vcf_path}...")
    commands = f"""
{TABIX} -f -p vcf {vcf_path}
"""
    script_name = f"{vcf_path}.index.slurm"
    job_name = f"index_vcf"

    create_slurm_script(script_name, job_name, commands, os.path.dirname(vcf_path))
    job_id = submit_to_slurm(script_name)
    print(f"Submitted VCF indexing job with Job ID: {job_id}")
    return job_id


def run_basevar(fq):
    for chromosome in PARAMETERS["chrs"]:
        # Step 1: Run BaseVar for the chromosome
        basevar_jobs = run_basevar_step(fq, chromosome)
        for job_id in basevar_jobs:
            wait_for_job_completion(job_id)

        # Step 2: Merge VCF files for the chromosome
        merge_job_id, merged_vcf = merge_vcf_files_slurm(fq, chromosome)
        wait_for_job_completion(merge_job_id)

        # Step 3: Index the merged VCF file
        index_job_id = index_vcf_file_slurm(merged_vcf)
        wait_for_job_completion(index_job_id)

        print(f"=== Completed processing for chromosome {chromosome} ===")

    print(f"=== Completed BaseVar pipeline for {fq} ===")
