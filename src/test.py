from helper.converter import convert_cram_to_fastq
from helper.path_define import cram_path, fastq_path
from helper.metrics import compare_fastq_sequences
import sys


def main():

    file1 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/0.1x/HG02019/sample_1/HG02019.fastq.gz"
    file2 = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/result/0.1x/HG02019/sample_2/HG02019.fastq.gz"
    compare_fastq_sequences(file1, file2)


if __name__ == "__main__":
    main()


