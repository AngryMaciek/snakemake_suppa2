# SUPPA2
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[General information about the project]

## Snakemake pipeline execution
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires Python 3 and can be most easily installed via the bioconda package from the anaconda cloud service.

### Step 1: Download and installation of Miniconda3
Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  source .bashrc
  ```

macOS:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  source .bashrc
  ```

### Step 2: Snakemake installation

Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest version:
  ```bash
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml 

Write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are four scripts to start the pipeline, depending on whether you want to run locally/on a SLURM computational cluster and if you prefer to use conda virtual environment or a specific container from dockerhub. In order to use singularity containers you need to have `singularity` installed. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  bash snakemake_local_run_singularity_container.sh
  bash snakemake_cluster_run_singularity_container.sh
  ```

## License

GPL v3.0


------------------

#############################################################
###############                               ###############
###############   SUPPA2 SNAKEMAKE PIPELINE   ###############
###############                               ###############
###############   Maciej Bak                  ###############
###############   maciej.bak@unibas.ch        ###############
###############                               ###############
#############################################################


1) Install miniconda and snakemake (Linux)

Snakemake is a workflow management system that helps to create and execute data processing pipelines.
It requires python3 and can be most easily installed via the bioconda package of the python anaconda distribution.
You're setup right after the following steps:

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source .bashrc

conda search -c bioconda snakemake
conda install -c bioconda snakemake=5.2.2

One will need pandas library as well, in case it's not already installed:

conda install -c anaconda pandas=0.23.4


2) Run the Pipeline

Before each run adjust the config.yaml with the right information (input/output/parameters).

* The workflow is rather simple and does not need a design table, all the information is in the config.
* The main input is a tsv file with transcript expression quantification (merged for all samples).
* IMPORTANT: The quantification has to be done with the same gtf file as specified for this pipeline.
* IMPORTANT: Sample names (columns of the main input) must match information specified in the config file.
* IMPORTANT: SUPPA2 differential splicing works only if there is >1 replicate in each condition.
* Differential Splicing is calculated as EXPERIMENT - CONTROL

There are three commands for running the pipeline:

* Write a DAG (directed acyclic graph) into dag.pdf and dag.png:
bash snakemake.dag_run

* Run the pipeline locally:
bash snakemake.local_run

* Run the pipeline on the computational cluster (SLURM workload menager):
(if required, adapt the 'slurm_cluster.json' and 'jobscript.sh' files before starting the run)
bash snakemake.cluster_run