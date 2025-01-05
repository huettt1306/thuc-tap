import os
from helper.config import PATHS
from helper.file_utils import extract_vcf

def base_dir(fq):
    return os.path.dirname(fq)

def samid(fq):
    return os.path.basename(fq).replace(".fastq.gz", "")

def tmp_outdir(fq):
    return os.path.join(base_dir(fq), "1tmp_files")

def batch1_final_outdir(fq):
    return os.path.join(base_dir(fq), "batch1_final_files")

def bamlist_dir(fq):
    return os.path.join(batch1_final_outdir(fq), f"{samid(fq)}.list")

def basevar_outdir(fq):
    return os.path.join(base_dir(fq), "basevar_output")

def basevar_vcf(fq, chromosome):
    return os.path.join(basevar_outdir(fq), f"NIPT_basevar_{chromosome}.vcf.gz")

def vcf_list_path(fq, chromosome):
    return os.path.join(basevar_outdir(fq), f"NIPT_basevar_{chromosome}.vcf.list")

def glimpse_outdir(fq):
    return os.path.join(base_dir(fq), "glimpse_output")

def vcf_prefix(chromosome):
    return f"CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}"

def get_vcf_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.filtered.shapeit2-duohmm-phased.vcf.gz")

def filtered_vcf_path(chromosome): 
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.biallelic.snp.maf0.001.sites.vcf.gz")

def filtered_tsv_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.biallelic.snp.maf0.001.sites.tsv.gz")

def chunks_path(chromosome):
    return os.path.join(PATHS["reference_path"], f"{vcf_prefix(chromosome)}.chunks.txt")

def glimpse_vcf(fq, chromosome):
    return os.path.join(glimpse_outdir(fq), "imputed_file_merged", f"glimpse.{chromosome}_imputed.vcf.gz")

def ground_truth_vcf(name, chromosome):
    path = os.path.join(PATHS["vcf_directory"], f"{name}_{chromosome}.vcf.gz")
    if os.path.exists(path):
        print(f"Ground truth VCF already exists: {path}")
        return path
    
    extract_vcf(name, path, chr=chromosome)
    return path

def statistic_outdir(fq, chromosome):
    return os.path.join(base_dir(fq), "statistic_output", f"{chromosome}")

def statistic_nipt_outdir(fq, chromosome, compare_with):
    return os.path.join(base_dir(fq), f"statistic_output_{compare_with}", f"{chromosome}")
