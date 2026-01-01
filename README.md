# UoM_CUT_and_RUN
Modularised code for CUT&amp;RUN pre-processing, designed specifically for University of Manchester CSF cluster


## How to use this pipeline
This pipeline is designed to take files from a raw fastq state, through various processing steps (see below), to aligned reads and called peaks.

1. File setup
As a general rule of thumb, raw files should be stored in a seperate directory to prevent from accidental editing of raw data, overwriting, or other administrative error. You can then generate symbolic links to raw data in your analysis directory.

Sample raw data should be contained within their own sample directory. For example:
- Data
-     Raw
-       Sample_1
-         Sample_1_Fastq_R1.fq.gz
-         Sample_1_Fastq_R2.fq.gz
-       Sample_2
-         Sample_2_Fastq_R1.fq.gz
-         Sample_2_Fastq_R2.fq.gz

2. Configuration file
This file contains pathways and parameters that will be largely specifc to invidiual analyses, and should be edited on a analysis-by-analysis basis. The main pre-requisite here is having access to the 'oa_functional_genomics" project directory. Some files a written to personal scratch to prevent unessessary storage usage.

3. Software requriments
Ideally, users should use the CUT_RUN conda environment. This is part of the "project system" implemented by Paul Martin within DMDS. Please contact jake.butler@manchester.ac.uk in the first instance, or Paul.Martin-2@manchester.ac.uk for details and setup. Software is then activated within SLURM job submission script via the following command:
- activate_project /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN.

*Please note i have also manually called software (via "module load") as i was pushed for time and didn't want to wait for mamba to install certain large software.*

Should further software be required, please activate this project from command line, and use mamba to add nessessary software.

4. SLURM job script
The job script is set up to run arrays (one array per sample). I make a sample.txt that contains one sample per line, this is then used to identify the "file to process" from the main.py and configuration.py.
I tend to run the pipeline in steps using the optional -s command as part of the pipeline. If ommitted, the entire pipeline runs.

# Input
Raw fastq files

# Output
