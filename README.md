
# swgs-binsize pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Philip Smith (@phil9s)

## Description

Generate absolute copy number profiles and signatures for multiple bin sizes from down sampled WGS data in order to compare feature distributions and quantitatively determine comparable segment sizes distributions between sWGS and SNP6.0 array copy number profiles.

Modified QDNAseq source code can be found [here](https://github.com/markowetzlab/QDNAseqmod)

## Table of contents

* [Pipeline setup](#pipeline-setup)
  + [Step 1 Clone the repo](#step-1-clone-the-repo)
  + [Step 2 Install conda](#step-2-install-conda)
    - [With Conda already installed](#with-conda-already-installed)
  + [Step 3 Installing additional dependencies](#step-3-installing-additional-dependencies)
  + [Step 4 Preparing the input files](#step-4-preparing-the-input-files)
    - [Sample sheet](#sample-sheet)
    - [config.yaml](#configyaml)
    - [profile configs](#profile-configs)
* [Running the pipeline](#running-the-pipeline)
* [Output](#output)

## Pipeline setup

### Step 1 Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
git clone https://github.com/Phil9S/swgs-binsize.git
cd swgs-binsize/
```

### Step 2 Install conda

Run the following to install conda whilst following the on-screen instructions.
- When asked to run `conda init` and initialise conda please respond with 'yes'

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $HOME/miniconda/
source ~/.bashrc
rm Miniconda3-latest-Linux-x86_64.sh
```
See [installing conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for more information.

#### With conda already installed

For systems where conda is already available the following requirements need to be met:
- conda must be available on the PATH
- conda version `4.8.3' or greater
- the location of the installation folder is required

Check the installed version of conda using the following:
```
conda -V
```
*If this command does not work then conda is also not available on the PATH*

Find your installation directory using the following:
```
whereis conda | sed 's%condabin/conda%%'
```

### Step 3 Installing additional dependencies

From within the repository directory, run the `install_env.sh` script to generate a conda environment and install custom packages:

```
./install_env.sh $HOME/miniconda/
```

If you used a previously installed conda build please use the conda or miniconda installation directory when running this section instead of '$HOME/miniconda/' to correctly initialise the conda environment.

The newly installed conda environment can be activated using the following:

```
conda activate swgs-binsize
```

### Step 4 Preparing the input files

#### Sample sheet

The workflow requires a single input file `sample_sheet.tsv` which is a tab-separated document detailing a number of ICGC/PCAWG identifiers. Edit this file to include the correct absolute path to the corresponding WGS tumour BAM file locations.

#### config.yaml

Edit the config.yaml (`config/config.yaml`) file to the desired output directory (include full directory path i.e. `/mnt/scratch/analysis/`).

#### profile configs

The profile configs (cluster_config.yaml & config.yaml) (`profile/Slurm/*`) contains the necessary information to configure the job submission parameters passed to [SLURM](https://slurm.schedmd.com/overview.html). This includes the number of concurrently summitted jobs, account/project name, partition/queue name, and default job resources (though these are low and should work on almost any cluster). Edit the values in these files as needed to suit the cluster environment and job limits.

## Running the pipeline

Once the pipeline and cluster parameters have been set and the `samplesheet.tsv` is prepared, the pipeline is ready to run.

With the `swgs-binsize` conda environment active, run the following:

Confirm the pipeline is configured correctly, run using the `dry-run` mode.

```
snakemake -n --profile profile/slurm/
```
If the previous step ran without error then run the following:
```
snakemake --profile profile/slurm/
```
This step may take some time to run depending on cluster availability and job run times, it is advised to run this in a `screen`.

## Output

The pipeline should generate a number of output folders containing absolute copy number profiles, copy number signature matrices, and lastly a folder containing the feature distributions.

The folder `sWGS_binsize` in the output directory contains the comparison outputs include distribution plots and statistical testing outputs. Segment size distributions are compared using a `mann whitney U` test, where a non-significant result demonstrates that segment size distribution of a given bin size is not significantly different from the segment size distribution generated by SNP6.0 array data. The bin size with a non-significant p-value is considered to be the best matching bin size to utilise for down sampling of PCAWG WGS samples.

## Licence
The contents of this repository are copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL). 

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

The methods implemented in the code are the subject of pending patent application GB 2114203.9.

Any commercial use of this code is prohibited.
