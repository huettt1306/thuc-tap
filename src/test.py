from helper.converter import convert_cram_to_fastq
from helper.path_define import cram_path, fastq_path
from helper.metrics import compare_fastq_sequences
import sys
from helper.file_utils import process_vcf

def main():
    file = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/vcf/HG02018_chr20.vcf.gz"
    process_vcf(file)


if __name__ == "__main__":
    main()


