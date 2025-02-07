import subprocess, os
from helper.config import TOOLS, PARAMETERS, PATHS
from helper.logger import setup_logger

# Thiết lập logger riêng cho quá trình chuyển đổi
logger = setup_logger(os.path.join(PATHS["logs"], "conversion.log"))

def convert_bam_to_fastq(bam_path, output_fastq_path):
    """
    Sử dụng samtools để chuyển đổi BAM sang FASTQ.GZ trực tiếp qua subprocess với số luồng tùy chỉnh.
    """
    logger.info(f"Converting BAM file {bam_path} to FASTQ.GZ at {output_fastq_path} using 8 threads...")

    # Lệnh samtools với số threads
    command = f"{TOOLS['samtools']} fastq -@ {PARAMETERS['threads']} -o {output_fastq_path} {bam_path}"

    try:
        # Gọi lệnh thông qua subprocess
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Kiểm tra kết quả
        if process.returncode != 0:
            error_message = f"Error converting BAM to FASTQ: {stderr.decode()}"
            logger.error(error_message)
            raise RuntimeError(error_message)
        else:
            logger.info(f"FASTQ file created successfully at {output_fastq_path}.")
    except Exception as e:
        logger.error(f"Failed to convert BAM to FASTQ: {e}")

def convert_cram_to_fastq(cram_path, output_fastq_path):
    """
    Sử dụng samtools để chuyển đổi CRAM sang FASTQ.GZ trực tiếp qua subprocess với số luồng tùy chỉnh.
    """
    logger.info(f"Converting CRAM file {cram_path} to FASTQ.GZ at {output_fastq_path} using 8 threads...")

    # Lệnh samtools với số threads
    command = f"{TOOLS['samtools']} fastq -@ {PARAMETERS['threads']} -o {output_fastq_path} {cram_path}"

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
            logger.info(f"FASTQ file created successfully at {output_fastq_path}.")
            os.remove(cram_path)
            logger.info(f"Original CRAM file {cram_path} deleted.")
    except Exception as e:
        logger.error(f"Failed to convert CRAM to FASTQ: {e}")


def convert_cram_to_fastq2(cram_path, output_fastq_gz_path):
    """
    Sử dụng samtools view để chuyển đổi CRAM sang FASTQ.GZ với pigz để tăng tốc nén.
    """
    logger.info(f"Converting CRAM file {cram_path} to FASTQ.GZ at {output_fastq_gz_path} using {PARAMETERS['threads']} threads...")

    # Lệnh samtools + pigz
    command = f"{TOOLS['samtools']} view -@ {PARAMETERS['threads']} -u {cram_path} | " \
              f"{TOOLS['samtools']} fastq - | {TOOLS['pigz']} -p {PARAMETERS['threads']} > {output_fastq_gz_path}"

    try:
        # Gọi subprocess để chạy lệnh
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Kiểm tra lỗi
        if process.returncode != 0:
            error_message = f"Error converting CRAM to FASTQ.GZ: {stderr.decode()}"
            logger.error(error_message)
            raise RuntimeError(error_message)
        else:
            logger.info(f"FASTQ.GZ file created successfully at {output_fastq_gz_path}.")
            os.remove(cram_path)  # Xóa file CRAM sau khi chuyển đổi
            logger.info(f"Original CRAM file {cram_path} deleted.")
    except Exception as e:
        logger.error(f"Failed to convert CRAM to FASTQ.GZ: {e}")
