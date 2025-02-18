from pipeline.generate import generate_single_sample, generate_nipt_sample
from pipeline.alignment import run_alignment_pipeline
from pipeline.basevar import run_basevar
from pipeline.glimpse import run_glimpse
from statistic.statistic import run_statistic

from pipeline.reference_panel_prepare import run_prepare_reference_panel
from helper.config import PARAMETERS, TRIO_DATA, PATHS
from helper.metrics import get_fastq_coverage
from helper.logger import setup_logger
from helper.file_utils import extract_lane1_fq
from helper.converter import convert_cram_to_fastq
from helper.path_define import fastq_path, fastq_path_lane1, fastq_path_lane2, cram_path, fastq_single_path, fastq_nipt_path
import os, sys
from concurrent.futures import ThreadPoolExecutor


logger = setup_logger(os.path.join(PATHS["logs"], "main.log"))


def pipeline_for_sample(fastq_dir):
    logger.info(f"Run all pipeline for sample in {fastq_dir}")
    run_alignment_pipeline(fastq_dir)
    run_basevar(fastq_dir)
    run_glimpse(fastq_dir)
    run_statistic(fastq_dir)

def prepare_data(name):
    print(f"Preparing data for {name}")
    if not os.path.exists(fastq_path_lane1(name)):
        convert_cram_to_fastq(cram_path(name), fastq_path_lane1(name), fastq_path_lane2(name))
    #return get_fastq_coverage(name)

def process_trio(trio_name, trio_info):
    """
    Xử lý một trio (bao gồm các bước pipeline cho từng mẫu).
    """
    logger.info(f"######## PROCESSING TRIO: {trio_name} ########")

    child_name = trio_info["child"]
    mother_name = trio_info["mother"]
    father_name = trio_info["father"]

    
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(prepare_data, [child_name, mother_name])

        #mother_avg_coverage = future_mother.result()
        #child_avg_coverage = future_child.result()


    for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
        logger.info(f"######## PROCESSING index {index} ########")

        for coverage in PARAMETERS["coverage"]:
            pipeline_for_sample(generate_single_sample(mother_name, coverage, index))

            for ff in PARAMETERS["ff"]:
                pipeline_for_sample(generate_nipt_sample(child_name, mother_name, father_name, coverage, ff, index))


def main():
    #run_prepare_reference_panel()
    if len(sys.argv) < 2:
        logger.error("Please provide a trio name to process.")
        sys.exit(1)
    
    trio_name = sys.argv[1]  
    if trio_name not in TRIO_DATA:
        logger.error(f"Trio {trio_name} not found in TRIO_DATA.")
        sys.exit(1)

    trio_info = TRIO_DATA[trio_name]
    process_trio(trio_name, trio_info)


if __name__ == "__main__":
    main()
