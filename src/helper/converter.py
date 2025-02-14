import subprocess, os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.logger import setup_logger

# Thiết lập logger riêng cho quá trình chuyển đổi
logger = setup_logger(os.path.join(PATHS["logs"], "conversion.log"))

def convert_cram_to_fastq(cram_path, output_fastq_path_1, output_fastq_path_2):
    """
    Sử dụng samtools để chuyển đổi CRAM sang FASTQ.GZ trực tiếp qua subprocess với số luồng tùy chỉnh.
    """
    logger.info(f"Converting CRAM file {cram_path} to FASTQ.GZ  using 8 threads...")

    # Lệnh samtools với số threads
    command = f"{TOOLS['samtools']} fastq -@ {PARAMETERS['threads']} -1 {output_fastq_path_1} -2 {output_fastq_path_2} {cram_path}"

    try:
        # Gọi lệnh thông qua subprocess
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Kiểm tra kết quả
        if process.returncode != 0:
            error_message = f"Error converting CRAM to FASTQ: {stderr.decode()}"
            logger.error(error_message)
            raise RuntimeError(error_message)
        else:
            logger.info(f"FASTQ file created successfully.")
            os.remove(cram_path)
            logger.info(f"Original CRAM file {cram_path} deleted.")
    except Exception as e:
        logger.error(f"Failed to convert CRAM to FASTQ: {e}")


def convert_genotype(genotype):
    if len(genotype) != 2:
        return "-1/-1"
    allen1 = genotype[0] if isinstance(genotype[0], int) else -1
    allen2 = genotype[1] if isinstance(genotype[1], int) else -1
    if allen1 > allen2:
        allen1, allen2 = allen2, allen1
    return f"{allen1}/{allen2}"
