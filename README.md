# Supplementary Scripts & Data

## Overview

This repository is an extension of the original **`bmds/CRISPR-benchmark`** project. It inherits the directory structure and adds new scripts and datasets to support the experiments presented in our paper.

*Note:* To ensure integrity and comparability with the original version, some files from `bmds/CRISPR-benchmark` are retained even if they are not directly involved in the computation pipeline.

## Prerequisites

The system requires **Miniconda3** to manage library packages and execution environments.

Due to differences in Python versions between tools (specifically the Azimuth model, which uses Python 2.x), the experimental workflow requires the setup of two independent virtual environments:

1. The `base` environment (updated original environment).
2. The `env2_conda` environment (new environment specifically for Azimuth).

## Installation & Environment Setup

Please execute the following commands to set up the execution environments:

```bash
# 1. Update the current base environment
conda env update --name base --file base_enviroment.yml

# 2. Create new env2_conda environment
conda env create -f env2_enviroment.yml

```

## Usage

Our main experimental scripts are designed to operate in the `base` environment. However, interacting with the Azimuth tool (based on Rule Set 2) requires the `env2_conda` environment.

* **Running main analysis scripts:**
```bash
conda activate base
# Execute pipelines here

```


* **Reproducing BMDS Benchmark results:**
To re-run the benchmark evaluation scripts from BMDS (limited to the results presented in our paper) on the new output data, please switch to the `env2_conda` environment:
```bash
conda activate env2_conda
# Execute BMDS scripts here

```



## Script Descriptions

Below is a list and the core functions of the scripts added to this repository:

* **`check_gw.py`**: Support script for extracting sgRNAs at corresponding regions from .csv files (output of `pipeline_GW.py`).
* **`experiment_functions_500.py`**: Contains core shared functions for calculation and data processing in experimental pipelines.
* **`get_rs2.py`**: Utility function to connect and retrieve data from the Azimuth model.
* **`pipeline_Doench.py`**: Automated pipeline for running accuracy evaluation experiments on the Doench et al. benchmark dataset.
* **`pipeline_Xu2015.py`**: Automated pipeline for running performance evaluation experiments on the selected dataset from Xu et al. (2015).
* **`pipeline_GW.py`**: Pipeline for whole-genome sgRNA screening experiments in Ecoli K-12.
* **`pipeline_cal_time.py`**: Pipeline for executing experiments and measuring time on the `500k` dataset. Output data is exported in `.normalised` format, fully standardized as input for BMDS evaluation scripts.
