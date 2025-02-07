from helper.converter import convert_cram_to_fastq, convert_cram_to_fastq2
from helper.path_define import cram_path, fastq_path
from helper.metrics import calculate_average_coverage
import sys


# Hàm để chuyển đổi một tệp
def convert(name):
    calculate_average_coverage(name)
    convert_cram_to_fastq(cram_path(name), fastq_path(name))
    
def main():

    if len(sys.argv) < 2:
        print("Please provide a name to process.")
        sys.exit(1)
    
    name = sys.argv[1]  # Lấy tên trio từ đầu vào dòng lệnh
    convert(name)


if __name__ == "__main__":
    main()


