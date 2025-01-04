import pysam
import os
from helper.file_utils import read_fastq_file


def calculate_average_coverage(name):
    """
    Tính độ phủ trung bình của dữ liệu reads.
    :param name: Tên mẫu (để hiển thị).
    """
    headers, reads, plus, qualities = read_fastq_file(name)
    genome_size = 3_200_000_000  # Kích thước genome tham chiếu (ví dụ: GRCh38)
    total_bases = sum(len(read) for read in reads)
    average_coverage = total_bases / genome_size

    print(f"Độ phủ trung bình của mẫu {name}: {average_coverage:.4f}")
    return average_coverage

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