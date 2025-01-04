import os
import random
from helper.file_utils import save_to_fastq, read_fastq_file
from helper.config import PATHS

generated_files = []

def select_reads_with_indices(data, qualities, headers, plus_separators, indices):
    """
    Chọn dữ liệu từ các mảng dựa trên chỉ mục ngẫu nhiên.
    """
    selected_reads = [data[i] for i in indices]
    selected_qualities = [qualities[i] for i in indices]
    selected_headers = [headers[i] for i in indices]
    selected_plus_separators = [plus_separators[i] for i in indices]
    return selected_reads, selected_qualities, selected_headers, selected_plus_separators

def generate_random_reads_files(data, qualities, headers, plus_separators, coverage, avg_coverage, output_prefix):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping creation.")
        return output_file

    ratio = coverage / avg_coverage
    # Tạo mảng chỉ mục và không xáo trộn
    indices = list(range(len(data)))

    # Chọn dữ liệu dựa trên tỷ lệ (độ dài dữ liệu * tỷ lệ)
    selected_indices = indices[:int(len(indices) * ratio)]

    selected_reads, selected_qualities, selected_headers, selected_plus_separators = select_reads_with_indices(
        data, qualities, headers, plus_separators, selected_indices
    )
    return save_to_fastq(output_file, selected_reads, selected_qualities, selected_headers, selected_plus_separators)

def generate_merge_files(
    child_data, child_qualities, child_headers, child_plus_separators,
    mother_data, mother_qualities, mother_headers, mother_plus_separators,
    coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix
):
    output_file = f"{output_prefix}.fastq.gz"
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping creation.")
        return output_file

    child_ratio = ff * coverage / child_avg_coverage
    mother_ratio = (1 - ff) * coverage / mother_avg_coverage

    # Tạo mảng chỉ mục và không xáo trộn
    child_indices = list(range(len(child_data)))
    mother_indices = list(range(len(mother_data)))

    # Chọn dữ liệu cho con và mẹ mà không xáo trộn
    child_selected_indices = child_indices[:int(len(child_indices) * child_ratio)]
    mother_selected_indices = mother_indices[:int(len(mother_indices) * mother_ratio)]

    # Lấy dữ liệu từ các chỉ mục đã chọn
    child_reads, child_qualities_selected, child_headers_selected, child_plus_selected = select_reads_with_indices(
        child_data, child_qualities, child_headers, child_plus_separators, child_selected_indices
    )
    mother_reads, mother_qualities_selected, mother_headers_selected, mother_plus_selected = select_reads_with_indices(
        mother_data, mother_qualities, mother_headers, mother_plus_separators, mother_selected_indices
    )

    combined_reads = child_reads + mother_reads
    combined_qualities = child_qualities_selected + mother_qualities_selected
    combined_headers = child_headers_selected + mother_headers_selected
    combined_plus = child_plus_selected + mother_plus_selected

    # Không cần xáo trộn mảng chỉ mục
    indices = list(range(len(combined_reads)))
    random.shuffle(indices)

    shuffled_reads = [combined_reads[i] for i in indices]
    shuffled_qualities = [combined_qualities[i] for i in indices]
    shuffled_headers = [combined_headers[i] for i in indices]
    shuffled_plus = [combined_plus[i] for i in indices]

    return save_to_fastq(output_file, shuffled_reads, shuffled_qualities, shuffled_headers, shuffled_plus)


def generate_single_sample(name, avg_coverage, coverage, index):
    """
    Tạo file dữ liệu tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    headers, reads, plus, qualities = read_fastq_file(name)
    sample_output_dir = os.path.join(PATHS["fastq_directory"], f"{coverage}x", name, f"sample_{index}")
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, name)
    output_dir = generate_random_reads_files(reads, qualities, headers, plus, coverage, avg_coverage, output_prefix)
    generated_files.append(output_dir)
    return output_dir


def generate_nipt_sample(child_name, child_avg_coverage, mother_name, mother_avg_coverage, coverage, ff, index):
    """
    Tạo file nipt tương ứng với cov và ind
    Return đường dẫn đến file fastq.gz tạo được
    """
    sample_output_dir = os.path.join(PATHS["fastq_directory"], f"{coverage}x", f"{child_name}_{mother_name}", f"{ff:.3f}", f"sample_{index}")
    os.makedirs(sample_output_dir, exist_ok=True)
    output_prefix = os.path.join(sample_output_dir, f"{child_name}_{mother_name}")
    output_dir = generate_merge_files(
        read_fastq_file(child_name), read_fastq_file(mother_name),
        coverage, child_avg_coverage, mother_avg_coverage, ff, output_prefix
    )
    generated_files.append(output_dir)
    return output_dir
 