import pysam, os
import subprocess
from helper.path_define import cram_path
from helper.config import TOOLS, PATHS

COVERAGE_FILE = os.path.join(PATHS["cram_directory"], "coverage.txt")

def calculate_average_coverage(name):
    # Đường dẫn đến file CRAM
    cram_file = cram_path(name)

    # Kiểm tra nếu kết quả đã tồn tại trong coverage.txt
    if os.path.exists(COVERAGE_FILE):
        with open(COVERAGE_FILE, 'r') as f:
            for line in f:
                sample, coverage = line.strip().split('\t')
                if sample == name:
                    # Nếu đã có kết quả, trả về ngay mà không tính lại
                    return float(coverage)

    # Tính toán độ bao phủ bằng lệnh samtools và awk
    try:
        cmd = f"{TOOLS['samtools']} depth {cram_file} | awk '{{total+=$3; count++}} END {{print total/count}}'"
        result = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        
        # Chuyển kết quả sang số thực
        average_coverage = float(result.stdout.strip())

        # Ghi kết quả vào file coverage.txt
        with open(COVERAGE_FILE, 'a') as f:
            f.write(f"{name}\t{average_coverage}\n")

        return average_coverage

    except subprocess.CalledProcessError as e:
        print(f"Error while running samtools or awk: {e.stderr}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None



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