from pipeline.generate import generate_single_sample, generate_nipt_sample
from pipeline.alignment import run_alignment_pipeline
from pipeline.basevar import run_basevar
from pipeline.glimpse import run_glimpse
from pipeline.statistic import run_statistic
from pipeline.reference_panel_prepare import prepare_reference_panel
from helper.config import PARAMETERS, TRIO_DATA
from helper.metrics import calculate_average_coverage


def pipeline_for_sample(fastq_dir):
    #run_alignment_pipeline(fastq_dir)
    #run_basevar(fastq_dir)
    run_glimpse(fastq_dir)
    run_statistic(fastq_dir)
    return


def main():
    #prepare_reference_panel()

    for trio_name, trio_info in TRIO_DATA.items():
        print(f"######## PROCESSING TRIO: {trio_name} ########")

        child_name = trio_info["child"]
        mother_name = trio_info["mother"]

        #child_avg_coverage = calculate_average_coverage(child_name)
        #mother_avg_coverage = calculate_average_coverage(mother_name)

        child_avg_coverage = 1
        mother_avg_coverage = 1

        for index in range(PARAMETERS["startSampleIndex"], PARAMETERS["endSampleIndex"] + 1):
            print(f"######## PROCESSING index {index} ########")

            for coverage in PARAMETERS["coverage"]:
                #pipeline_for_sample(generate_single_sample(child_name, child_avg_coverage, coverage, index))
                #pipeline_for_sample(generate_single_sample(mother_name, mother_avg_coverage, coverage, index))
                
                for ff in PARAMETERS["ff"]:
                    pipeline_for_sample(generate_nipt_sample(child_name, child_avg_coverage, mother_name, mother_avg_coverage, coverage, ff, index))            


if __name__ == "__main__":
    main()
