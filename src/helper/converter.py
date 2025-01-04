from helper.config import TOOLS
from helper.slurm import create_slurm_script, submit_to_slurm, wait_for_job_completion


def convert_bam_to_fastq(bam_path, output_fastq_path):
    """
    Sử dụng samtools để chuyển đổi BAM sang FASTQ.GZ.
    """
    print(f"Converting BAM file {bam_path} to FASTQ.GZ at {output_fastq_path}...")

    commands = f"""
{TOOLS['samtools']} fastq -o {output_fastq_path} {bam_path}
    """

    script_name = f"convert_bam_to_fastq_{bam_path.split('/')[-1].split('.')[0]}.slurm"
    job_name = f"convert_bam_to_fastq"

    create_slurm_script(script_name, job_name, commands, output_dir=output_fastq_path)
    job_id = submit_to_slurm(script_name)
    wait_for_job_completion(job_id)

    print(f"FASTQ file created at {output_fastq_path}.")

# TODO: đoạn mã convert bên dưới chưa đúng
def convert_cram_to_fastq(cram_path, output_fastq_path):
    """
    Sử dụng samtools để chuyển đổi CRAM sang FASTQ.GZ, sử dụng SLURM.
    """
    print(f"Converting CRAM file {cram_path} to FASTQ.GZ at {output_fastq_path}...")

    commands = f"""
{TOOLS['samtools']} fastq -o {output_fastq_path} {cram_path}
    """

    script_name = f"convert_cram_to_fastq_{cram_path.split('/')[-1].split('.')[0]}.slurm"
    job_name = f"convert_cram_to_fastq"

    create_slurm_script(script_name, job_name, commands, output_dir=output_fastq_path)
    job_id = submit_to_slurm(script_name)
    wait_for_job_completion(job_id)

    print(f"FASTQ file created at {output_fastq_path}.")