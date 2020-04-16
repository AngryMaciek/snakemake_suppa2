# Snakemake pipeline for Alternative Splicing analysis with SUPPA2
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[SUPPA2](https://github.com/comprna/SUPPA) is a very nice tool to quantify AS events from RNA-Seq samples and perform differential splicing analysis.
This repository is a very small snakemake workflow that I use for automated and reproducible analyses in my reseach.

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

### Step 2: Pandas and Snakemake installation

To execute the workflow one would require pandas python library and snakemake workflow menager.  
Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest versions:
  ```bash
  conda install -c conda-forge pandas
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml  
Please note that: 

* The main input is a tsv file with transcript expression quantification (merged for all samples). The quantification has to be performed with the same gtf file as specified for this pipeline.
* All sample names (columns of the main input) have to match information specified in the config file.
* Differential splicing analysis with SUPPA2 works only if there is >1 replicate in each condition.
* Differential Splicing is calculated as EXPERIMENT - CONTROL groups

Once the metadata are ready write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are two scripts to start the pipeline, depending on whether you want to run locally or on a SLURM computational cluster. In order to execute SUPPA2 snakemake automatically creates an internal conda virtual environment and installs the tool from anaconda cloud service. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  ```

## License

Apache 2.0
