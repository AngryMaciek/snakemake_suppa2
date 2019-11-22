##############################################################################
#
#   Snakemake pipeline:
#   SUPPA2
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 22-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################

# imports
import sys
import os

# local rules
localrules: create_output_dir, all

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TXT_final_results = \
            expand(os.path.join("{output_dir}", "results.txt"),
                output_dir=config["output_dir"])

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_random_samples = os.path.join("{output_dir}", "random_samples"),
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_random_samples}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Sample some random data
##############################################################################

rule generate_files:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        SCRIPT = \
            os.path.join(config["src_dir"], "mb_random_sample.py")
    output:
        TXT_random_sample = \
            os.path.join("{output_dir}", "random_samples", "{file}")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "generate_files_{file}.log"),
        queue = "30min",
        time = "0:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "generate_files_{file}.log"),
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "generate_files_{file}_benchmark.log")
    conda:
        "packages.yaml"
    singularity:
        ""
    shell:
        """
        python {input.SCRIPT} \
        --outfile {output.TXT_random_sample} \
        &> {log.LOG_local_log};
        """

##############################################################################
### Merge the results
##############################################################################

rule merge_results:
    input:
        TXT_result_files = \
            lambda wildcards: [os.path.join(wildcards.output_dir,
                "random_samples", f) for f in config["samples_filenames"]]
    output:
        TXT_final_results = os.path.join("{output_dir}", "results.txt")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/merge_results.log"),
        queue = "30min",
        time = "00:05:00"
    resources:
        threads = 1,
        mem = 5000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_results.log")
    run:
        # read all the sampled numbers:
        numbers = []
        for i in input.TXT_result_files:
            with open(i) as f:
                numbers.append(f.read().splitlines()[0])
        # save into one file:
        with open(output.TXT_final_results, "w") as outfile:
                outfile.write("\n".join(numbers))

