import pysam, os
import subprocess
from helper.path_define import fastq_path_land1
from helper.config import TOOLS, PATHS
from helper.logger import setup_logger


COVERAGE_FILE = os.path.join(PATHS["fastq_directory"], "coverage.txt")
logger = setup_logger(os.path.join(PATHS["logs"], "metrics.log"))

def get_fastq_coverage(name):
    """
    Tính coverage của file FASTQ
    """
    # Kiểm tra nếu kết quả đã tồn tại trong coverage.txt
    if os.path.exists(COVERAGE_FILE):
        with open(COVERAGE_FILE, 'r') as f:
            for line in f:
                sample, coverage = line.strip().split('\t')
                if sample == name:
                    # Nếu đã có kết quả, trả về ngay mà không tính lại
                    return float(coverage)


    # Chạy `seqkit stats` để lấy tổng số base (sum_len)
    input_fastq = fastq_path_land1(name)
    result = subprocess.run(
        f"{TOOLS['seqkit']} stats {input_fastq} | tail -n 1",
        shell=True,
        capture_output=True,
        text=True
    )

    # Phân tích kết quả đầu ra từ `seqkit stats`
    stats_values = result.stdout.split()
    
    # Kiểm tra xem `seqkit` có chạy thành công không
    if len(stats_values) < 5:
        raise ValueError("Không thể lấy thống kê từ seqkit. Hãy kiểm tra file FASTQ.")

    # Lấy tổng số nucleotide (sum_len) từ cột 5
    total_bases = int(stats_values[4].replace(",", ""))  

    # Tính coverage
    genome_size = 3200000000
    coverage = total_bases / genome_size

    with open(COVERAGE_FILE, 'a') as f:
        f.write(f"{name}\t{coverage}\n")

    return coverage


def evaluate_vcf(ground_truth_file, test_file):
    """
    So sánh hai file VCF và tính toán các chỉ số.
    """
    gt_vcf = pysam.VariantFile(ground_truth_file)
    test_vcf = pysam.VariantFile(test_file)

    ground_truth_variants = {(rec.chrom, rec.pos, rec.ref, tuple(rec.alts)) for rec in gt_vcf.fetch()}
    test_variants = {(rec.chrom, rec.pos, rec.ref, tuple(rec.alts)) for rec in test_vcf.fetch()}

    true_positives = test_variants & ground_truth_variants
    false_positives = test_variants - ground_truth_variants
    false_negatives = ground_truth_variants - test_variants

    total_ground_truth = len(ground_truth_variants)
    total_test_variants = len(test_variants)

    accuracy = len(true_positives) / total_test_variants * 100 if total_test_variants > 0 else 0
    precision = len(true_positives) / (len(true_positives) + len(false_positives)) * 100 if (len(true_positives) + len(false_positives)) > 0 else 0
    recall = len(true_positives) / total_ground_truth * 100 if total_ground_truth > 0 else 0

    return {
        "total_ground_truth": total_ground_truth,
        "total_test_variants": total_test_variants,
        "true_positives": len(true_positives),
        "false_positives": len(false_positives),
        "false_negatives": len(false_negatives),
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
    }