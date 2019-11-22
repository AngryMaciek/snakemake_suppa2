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
import pandas as pd

# local rules
localrules: all, create_result_dir, extract_expression, \
    merge_SE, merge_MX, merge_RI, merge_AF, merge_AL, merge_A3, merge_A5, \
    merge_isoform, merge_differential_splicing

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TSV_SE_ALL = expand(os.path.join("{output_dir}", "SE_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_RI_ALL = expand(os.path.join("{output_dir}", "RI_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_MX_ALL = expand(os.path.join("{output_dir}", "MX_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_AF_ALL = expand(os.path.join("{output_dir}", "AF_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_AL_ALL = expand(os.path.join("{output_dir}", "AL_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_A3_ALL = expand(os.path.join("{output_dir}", "A3_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_A5_ALL = expand(os.path.join("{output_dir}", "A5_ALL.psi"), \
            output_dir=config["result_dir"]),
        TSV_isoform_psivec = \
            expand(os.path.join("{output_dir}", "Delta_psi_isoform.psivec"), \
                output_dir=config["result_dir"]),
        TSV_isoform_dpsi = \
            expand(os.path.join("{output_dir}", "Delta_psi_isoform.dpsi"), \
                output_dir=config["result_dir"]),
        TSV_merged_AS_table = \
            expand(os.path.join("{output_dir}", "local_AS.tsv"), \
                output_dir=config["result_dir"])

##############################################################################
### Create directories for the result
##############################################################################

rule create_result_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_main_dir = "{output_dir}",
        LOG_cluster_log = "{output_dir}/cluster_log",
    log:
        DIR_local_log = "{output_dir}/local_log"
    shell:
        """
        mkdir -p {params.DIR_main_dir}; \
        mkdir -p {params.LOG_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### SUPPA: prepare local AS event type and transcript events tables
##############################################################################

rule suppa_generate_local_events:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created")
    output:
        TSV_SE_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_SE_variable_"+config["threshold"]+".ioe"),
        TSV_A5_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A5_variable_"+config["threshold"]+".ioe"),
        TSV_A3_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A3_variable_"+config["threshold"]+".ioe"),
        TSV_RI_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_RI_variable_"+config["threshold"]+".ioe"),
        TSV_MX_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_MX_variable_"+config["threshold"]+".ioe"),
        TSV_AF_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AF_variable_"+config["threshold"]+".ioe"),
        TSV_AL_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AL_variable_"+config["threshold"]+".ioe")
    params:
        GTF_gtf = config["gtf"],
        threshold = config["threshold"],
        # workaround due to suppa file creation:
        STRING_output_prefix = "{output_dir}/AS_events/events",
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_generate_local_events.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_generate_local_events.log"),
    resources:
        threads = 1,
        mem = 10000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_generate_local_events_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py generateEvents \
        --input-file {params.GTF_gtf} \
        --format ioe \
        --output-file {params.STRING_output_prefix} \
        --event-type {{SE,SS,MX,RI,FL}} \
        --boundary V \
        --threshold {params.threshold} \
        &> {log.LOG_local_log}
        """

rule suppa_generate_transcript_events:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created")
    output:
        TSV_transcript_events_ioi = \
            os.path.join("{output_dir}", "AS_events", "transcript_events.ioi")
    params:
        GTF_gtf = config["gtf"],
        # workaround due to suppa file creation:
        STRING_output_prefix = \
            os.path.join("{output_dir}", "AS_events", "transcript_events"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_generate_transcript_events.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_generate_transcript_events.log"),
    resources:
        threads = 1,
        mem = 10000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_generate_transcript_events_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py generateEvents \
        --input-file {params.GTF_gtf} \
        --format ioi \
        --output-file {params.STRING_output_prefix} \
        &> {log.LOG_local_log}
        """

##############################################################################
### Extract transcripts expression per sample
##############################################################################

rule extract_expression:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created")
    output:
        DIR_quantification_dir = \
            directory(os.path.join("{output_dir}", "quantification"))
    params:
        TSV_transcript_quantification = config["transcript_quantification"],
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "extract_expression.log"),
    run:
        # split transcript expression table column-wise and condition-wise
        os.mkdir(output.DIR_quantification_dir)
        quant = pd.read_csv(config["transcript_quantification"], \
            sep="\t", index_col=0)
        for s in quant.columns.values:
            outfile_name = \
                os.path.join(output.DIR_quantification_dir, s+".tsv")
            quant[[s]].to_csv(outfile_name, sep="\t", index_label=False)
        temp_str = \
            os.path.join(output.DIR_quantification_dir,"CONTROL.tsv")
        quant[params.LIST_control_samples].to_csv(temp_str, \
            sep="\t", index_label=False)
        temp_str = \
            os.path.join(output.DIR_quantification_dir,"EXPERIMENT.tsv")
        quant[params.LIST_experiment_samples].to_csv(temp_str, \
            sep="\t", index_label=False)

##############################################################################
### SUPPA: PSI calculation per local events
##############################################################################

rule suppa_psi_per_event_SE:
    input:
        TSV_SE_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_SE_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "SE.psi")
    params:
        STRING_output_prefix = os.path.join("{output_dir}", "{sample}", "SE"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_SE_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_SE_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_SE_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.SE_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_MX:
    input:
        TSV_MX_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_MX_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "MX.psi")
    params:
        STRING_output_prefix = os.path.join("{output_dir}", "{sample}", "MX"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_MX_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_MX_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_MX_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.MX_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_RI:
    input:
        TSV_RI_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_RI_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "RI.psi")
    params:
        STRING_output_prefix = os.path.join("{output_dir}", "{sample}", "RI"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_RI_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_RI_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_RI_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.RI_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_AF:
    input:
        TSV_AF_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AF_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "AF.psi")
    params:
        STAFNG_output_prefix = os.path.join("{output_dir}", "{sample}", "AF"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_AF_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_AF_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_AF_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.AF_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_AL:
    input:
        TSV_AL_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AL_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "AL.psi")
    params:
        STALNG_output_prefix = os.path.join("{output_dir}", "{sample}", "AL"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_AL_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_AL_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_AL_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.AL_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_A3:
    input:
        TSV_A3_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A3_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "A3.psi")
    params:
        STA3NG_output_prefix = os.path.join("{output_dir}", "{sample}", "A3"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_A3_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_A3_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_A3_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.A3_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

rule suppa_psi_per_event_A5:
    input:
        TSV_A5_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A5_variable_"+config["threshold"]+".ioe"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = os.path.join("{output_dir}", "{sample}", "A5.psi")
    params:
        STA5NG_output_prefix = os.path.join("{output_dir}", "{sample}", "A5"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        tpm_filter = config["total_filter"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_event_A5_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_event_A5_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_event_A5_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerEvent \
        --ioe-file {input.A5_ioe} \
        --expression-file {params.quant} \
        --total-filter {params.tpm_filter} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log}
        """

##############################################################################
### SUPPA: PSI calculation per transcript events
##############################################################################

rule suppa_psi_per_isoform:
    input:
        TSV_transcript_events_ioi = \
            os.path.join("{output_dir}", "AS_events", \
                "transcript_events.ioi"),
        DIR_quantification_dir = \
            os.path.join("{output_dir}", "quantification")
    output:
        TSV_PSI_output = \
            os.path.join("{output_dir}", "{sample}", "_isoform.psi")
    params:
        GTF_gtf = config["gtf"],
        STRING_output_prefix = os.path.join("{output_dir}", "{sample}"),
        TSV_sample_quant = \
            os.path.join("{output_dir}", "quantification", "{sample}.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_psi_per_isoform_{sample}.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "suppa_psi_per_isoform_{sample}.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_psi_per_isoform_{sample}_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py psiPerIsoform \
        --gtf {params.GTF_gtf} \
        --expression-file {params.TSV_sample_quant} \
        --output-file {params.STRING_output_prefix} \
        &> {log.LOG_local_log}
        """

##############################################################################
### Merge multiple tables
##############################################################################

rule merge_SE:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "SE.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_SE_CONTROL = os.path.join("{output_dir}", "SE_CONTROL.psi"),
        TSV_SE_EXPERIMENT = os.path.join("{output_dir}", "SE_EXPERIMENT.psi"),
        TSV_SE_ALL = os.path.join("{output_dir}", "SE_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_SE.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"SE.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"SE_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"SE.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"SE_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"SE.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"SE_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_MX:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "MX.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_MX_CONTROL = os.path.join("{output_dir}", "MX_CONTROL.psi"),
        TSV_MX_EXPERIMENT = os.path.join("{output_dir}", "MX_EXPERIMENT.psi"),
        TSV_MX_ALL = os.path.join("{output_dir}", "MX_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_MX.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"MX.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"MX_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"MX.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"MX_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"MX.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"MX_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_RI:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "RI.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_RI_CONTROL = os.path.join("{output_dir}", "RI_CONTROL.psi"),
        TSV_RI_EXPERIMENT = os.path.join("{output_dir}", "RI_EXPERIMENT.psi"),
        TSV_RI_ALL = os.path.join("{output_dir}", "RI_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_RI.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"RI.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"RI_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"RI.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"RI_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"RI.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"RI_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_AF:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "AF.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_AF_CONTROL = os.path.join("{output_dir}", "AF_CONTROL.psi"),
        TSV_AF_EXPERIMENT = os.path.join("{output_dir}", "AF_EXPERIMENT.psi"),
        TSV_AF_ALL = os.path.join("{output_dir}", "AF_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_AF.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"AF.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AF_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"AF.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AF_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"AF.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AF_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_AL:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "AL.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_AL_CONTROL = os.path.join("{output_dir}", "AL_CONTROL.psi"),
        TSV_AL_EXPERIMENT = os.path.join("{output_dir}", "AL_EXPERIMENT.psi"),
        TSV_AL_ALL = os.path.join("{output_dir}", "AL_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_AL.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"AL.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AL_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"AL.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AL_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"AL.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"AL_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_A3:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "A3.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_A3_CONTROL = os.path.join("{output_dir}", "A3_CONTROL.psi"),
        TSV_A3_EXPERIMENT = os.path.join("{output_dir}", "A3_EXPERIMENT.psi"),
        TSV_A3_ALL = os.path.join("{output_dir}", "A3_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_A3.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"A3.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A3_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"A3.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A3_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"A3.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A3_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_A5:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "A5.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_A5_CONTROL = os.path.join("{output_dir}", "A5_CONTROL.psi"),
        TSV_A5_EXPERIMENT = os.path.join("{output_dir}", "A5_EXPERIMENT.psi"),
        TSV_A5_ALL = os.path.join("{output_dir}", "A5_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_A5.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"A5.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A5_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"A5.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A5_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"A5.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"A5_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

rule merge_isoform:
    input:
        expand(os.path.join("{output_dir}", "{sample}", "_isoform.psi"), \
            output_dir=config["result_dir"], \
            sample=config["CONTROL"]+config["EXPERIMENT"])
    output:
        TSV_isoform_CONTROL = \
            os.path.join("{output_dir}", "isoform_CONTROL.psi"),
        TSV_isoform_EXPERIMENT = \
            os.path.join("{output_dir}", "isoform_EXPERIMENT.psi"),
        TSV_isoform_ALL = \
            os.path.join("{output_dir}", "isoform_ALL.psi"),
    params:
        DIR_output_dir = "{output_dir}",
        LIST_control_samples = config["CONTROL"],
        LIST_experiment_samples = config["EXPERIMENT"],
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_isoform.log")
    run:
        # merge PSI tables for both conditions and for all samples
        df_list = []
        for i in params.control_samples:
            fname = os.path.join(params.output_dir,i,"_isoform.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"isoform_CONTROL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"_isoform.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"isoform_EXPERIMENT.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)
        df_list = []
        for i in params.control_samples+params.experiment_samples:
            fname = os.path.join(params.output_dir,i,"_isoform.psi")
            df_list.append(pd.read_csv(fname, sep="\t", index_col=0))
        merged_df = pd.concat(df_list, axis=1).fillna("nan")
        fname = os.path.join(params.output_dir,"isoform_ALL.psi")
        merged_df.to_csv(fname, sep="\t", index_label=False)

##############################################################################
### Differential splicing analysis for local AS events
##############################################################################

rule suppa_diffSplice_SE:
    input:
        TSV_SE_CONTROL = os.path.join("{output_dir}", "SE_CONTROL.psi"),
        TSV_SE_EXPERIMENT = os.path.join("{output_dir}", "SE_EXPERIMENT.psi"),
    output:
        TSV_SE_psivec = os.path.join("{output_dir}", "Delta_psi_SE.psivec"),
        TSV_SE_dpsi = os.path.join("{output_dir}", "Delta_psi_SE.dpsi")
    params:
        TSV_SE_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_SE_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STRING_output_prefix = "Delta_psi_SE",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_SE.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_SE.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_SE_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.SE_ioe} \
        --psi {input.SE_CONTROL} {input.SE_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_SE.psivec {params.output_dir}/Delta_psi_SE.psivec &&
        mv Delta_psi_SE.dpsi {params.output_dir}/Delta_psi_SE.dpsi;
        """

rule suppa_diffSplice_MX:
    input:
        TSV_MX_CONTROL = os.path.join("{output_dir}", "MX_CONTROL.psi"),
        TSV_MX_EXPERIMENT = os.path.join("{output_dir}", "MX_EXPERIMENT.psi"),
    output:
        TSV_MX_psivec = os.path.join("{output_dir}", "Delta_psi_MX.psivec"),
        TSV_MX_dpsi = os.path.join("{output_dir}", "Delta_psi_MX.dpsi")
    params:
        TSV_MX_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_MX_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STRING_output_prefix = "Delta_psi_MX",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_MX.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_MX.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_MX_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.MX_ioe} \
        --psi {input.MX_CONTROL} {input.MX_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_MX.psivec {params.output_dir}/Delta_psi_MX.psivec &&
        mv Delta_psi_MX.dpsi {params.output_dir}/Delta_psi_MX.dpsi;
        """

rule suppa_diffSplice_RI:
    input:
        TSV_RI_CONTROL = os.path.join("{output_dir}", "RI_CONTROL.psi"),
        TSV_RI_EXPERIMENT = os.path.join("{output_dir}", "RI_EXPERIMENT.psi"),
    output:
        TSV_RI_psivec = os.path.join("{output_dir}", "Delta_psi_RI.psivec"),
        TSV_RI_dpsi = os.path.join("{output_dir}", "Delta_psi_RI.dpsi")
    params:
        TSV_RI_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_RI_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STRING_output_prefix = "Delta_psi_RI",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_RI.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_RI.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_RI_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.RI_ioe} \
        --psi {input.RI_CONTROL} {input.RI_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_RI.psivec {params.output_dir}/Delta_psi_RI.psivec &&
        mv Delta_psi_RI.dpsi {params.output_dir}/Delta_psi_RI.dpsi;
        """

rule suppa_diffSplice_AF:
    input:
        TSV_AF_CONTROL = os.path.join("{output_dir}", "AF_CONTROL.psi"),
        TSV_AF_EXPERIMENT = os.path.join("{output_dir}", "AF_EXPERIMENT.psi"),
    output:
        TSV_AF_psivec = os.path.join("{output_dir}", "Delta_psi_AF.psivec"),
        TSV_AF_dpsi = os.path.join("{output_dir}", "Delta_psi_AF.dpsi")
    params:
        TSV_AF_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AF_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STAFNG_output_prefix = "Delta_psi_AF",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_AF.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_AF.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_AF_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.AF_ioe} \
        --psi {input.AF_CONTROL} {input.AF_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_AF.psivec {params.output_dir}/Delta_psi_AF.psivec &&
        mv Delta_psi_AF.dpsi {params.output_dir}/Delta_psi_AF.dpsi;
        """

rule suppa_diffSplice_AL:
    input:
        TSV_AL_CONTROL = os.path.join("{output_dir}", "AL_CONTROL.psi"),
        TSV_AL_EXPERIMENT = os.path.join("{output_dir}", "AL_EXPERIMENT.psi"),
    output:
        TSV_AL_psivec = os.path.join("{output_dir}", "Delta_psi_AL.psivec"),
        TSV_AL_dpsi = os.path.join("{output_dir}", "Delta_psi_AL.dpsi")
    params:
        TSV_AL_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_AL_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STALNG_output_prefix = "Delta_psi_AL",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_AL.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_AL.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_AL_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.AL_ioe} \
        --psi {input.AL_CONTROL} {input.AL_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_AL.psivec {params.output_dir}/Delta_psi_AL.psivec &&
        mv Delta_psi_AL.dpsi {params.output_dir}/Delta_psi_AL.dpsi;
        """

rule suppa_diffSplice_A3:
    input:
        TSV_A3_CONTROL = os.path.join("{output_dir}", "A3_CONTROL.psi"),
        TSV_A3_EXPERIMENT = os.path.join("{output_dir}", "A3_EXPERIMENT.psi"),
    output:
        TSV_A3_psivec = os.path.join("{output_dir}", "Delta_psi_A3.psivec"),
        TSV_A3_dpsi = os.path.join("{output_dir}", "Delta_psi_A3.dpsi")
    params:
        TSV_A3_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A3_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STA3NG_output_prefix = "Delta_psi_A3",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_A3.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_A3.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_A3_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.A3_ioe} \
        --psi {input.A3_CONTROL} {input.A3_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_A3.psivec {params.output_dir}/Delta_psi_A3.psivec &&
        mv Delta_psi_A3.dpsi {params.output_dir}/Delta_psi_A3.dpsi;
        """

rule suppa_diffSplice_A5:
    input:
        TSV_A5_CONTROL = os.path.join("{output_dir}", "A5_CONTROL.psi"),
        TSV_A5_EXPERIMENT = os.path.join("{output_dir}", "A5_EXPERIMENT.psi"),
    output:
        TSV_A5_psivec = os.path.join("{output_dir}", "Delta_psi_A5.psivec"),
        TSV_A5_dpsi = os.path.join("{output_dir}", "Delta_psi_A5.dpsi")
    params:
        TSV_A5_ioe = os.path.join("{output_dir}", "AS_events", \
            "events_A5_variable_"+config["threshold"]+".ioe"),
        DIR_output_dir = "{output_dir}",
        STA5NG_output_prefix = "Delta_psi_A5",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_A5.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_A5.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_A5_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.A5_ioe} \
        --psi {input.A5_CONTROL} {input.A5_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_A5.psivec {params.output_dir}/Delta_psi_A5.psivec &&
        mv Delta_psi_A5.dpsi {params.output_dir}/Delta_psi_A5.dpsi;
        """

##############################################################################
### Differential splicing analysis for transcript events
##############################################################################

rule suppa_diffSplice_transcripts:
    input:
        TSV_isoform_CONTROL = \
            os.path.join("{output_dir}", "isoform_CONTROL.psi"),
        TSV_isoform_EXPERIMENT = \
            os.path.join("{output_dir}", "isoform_EXPERIMENT.psi"),
    output:
        TSV_isoform_psivec = \
            os.path.join("{output_dir}", "Delta_psi_isoform.psivec"),
        TSV_isoform_dpsi = \
            os.path.join("{output_dir}", "Delta_psi_isoform.dpsi")
    params:
        TSV_transcript_events_ioi = \
            os.path.join("{output_dir}", "AS_events", \
                "transcript_events.ioi"),
        DIR_output_dir = "{output_dir}",
        method = "empirical",
        area = "1000",
        lower_bound = config["lower_bound"],
        gene_correction = "--gene-correction",
        STRING_output_prefix = "Delta_psi_isoform",
        TSV_control_quant = \
            os.path.join("{output_dir}", "quantification", "CONTROL.tsv"),
        TSV_experiment_quant = \
            os.path.join("{output_dir}", "quantification", "EXPERIMENT.tsv"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "suppa_diffSplice_transcripts.log"),
        queue = "6hours",
        time = "06:00:00"
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "suppa_diffSplice_transcripts.log")
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}", \
            "cluster_log", "suppa_diffSplice_transcripts_benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        suppa.py diffSplice \
        --method {params.method} \
        --input {params.ioi} \
        --psi {input.isoform_CONTROL} {input.isoform_EXPERIMENT} \
        --tpm {params.control_quant} {params.experiment_quant} \
        --area {params.area} \
        -gc \
        --lower-bound {params.lower_bound} \
        {params.gene_correction} \
        --output {params.output_prefix} \
        &> {log.LOG_local_log} &&
        sleep 5 &&
        mv Delta_psi_isoform.psivec \
        {params.output_dir}/Delta_psi_isoform.psivec &&
        mv Delta_psi_isoform.dpsi \
        {params.output_dir}/Delta_psi_isoform.dpsi;
        """

##############################################################################
### Merge differential splicing analysis results of local events
##############################################################################

rule merge_differential_splicing:
    input:
        TSV_SE_psivec = os.path.join("{output_dir}", "Delta_psi_SE.psivec"),
        TSV_SE_dpsi = os.path.join("{output_dir}", "Delta_psi_SE.dpsi"),
        TSV_MX_psivec = os.path.join("{output_dir}", "Delta_psi_MX.psivec"),
        TSV_MX_dpsi = os.path.join("{output_dir}", "Delta_psi_MX.dpsi"),
        TSV_RI_psivec = os.path.join("{output_dir}", "Delta_psi_RI.psivec"),
        TSV_RI_dpsi = os.path.join("{output_dir}", "Delta_psi_RI.dpsi"),
        TSV_AF_psivec = os.path.join("{output_dir}", "Delta_psi_AF.psivec"),
        TSV_AF_dpsi = os.path.join("{output_dir}", "Delta_psi_AF.dpsi"),
        TSV_AL_psivec = os.path.join("{output_dir}", "Delta_psi_AL.psivec"),
        TSV_AL_dpsi = os.path.join("{output_dir}", "Delta_psi_AL.dpsi"),
        TSV_A3_psivec = os.path.join("{output_dir}", "Delta_psi_A3.psivec"),
        TSV_A3_dpsi = os.path.join("{output_dir}", "Delta_psi_A3.dpsi"),
        TSV_A5_psivec = os.path.join("{output_dir}", "Delta_psi_A5.psivec"),
        TSV_A5_dpsi = os.path.join("{output_dir}", "Delta_psi_A5.dpsi"),
    output:
        TSV_merged_AS_table = os.path.join("{output_dir}", "local_AS.tsv")
    log:
        LOG_local_log = os.path.join("{output_dir}", "local_log", \
            "merge_differential_splicing.log")
    run:
        # merge all local event types:
        df_list = []
        SE_df = pd.read_csv(input.SE_dpsi, sep="\t")
        SE_df.columns = ["dPSI","pval"]
        df_list.append(SE_df)
        MX_df = pd.read_csv(input.MX_dpsi, sep="\t")
        MX_df.columns = ["dPSI","pval"]
        df_list.append(MX_df)
        RI_df = pd.read_csv(input.RI_dpsi, sep="\t")
        RI_df.columns = ["dPSI","pval"]
        df_list.append(RI_df)
        AF_df = pd.read_csv(input.AF_dpsi, sep="\t")
        AF_df.columns = ["dPSI","pval"]
        df_list.append(AF_df)
        AL_df = pd.read_csv(input.AL_dpsi, sep="\t")
        AL_df.columns = ["dPSI","pval"]
        df_list.append(AL_df)
        A5_df = pd.read_csv(input.A5_dpsi, sep="\t")
        A5_df.columns = ["dPSI","pval"]
        df_list.append(A5_df)
        A3_df = pd.read_csv(input.A3_dpsi, sep="\t")
        A3_df.columns = ["dPSI","pval"]
        df_list.append(A3_df)
        merged_dpsi = pd.concat(df_list)
        merged_dpsi.to_csv(output.merged_AS_table, sep="\t")

