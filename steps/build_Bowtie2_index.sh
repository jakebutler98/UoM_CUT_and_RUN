#!/bin/bash --login
#SBATCH -J index_ecoli                 # Job name
#SBATCH -p multicore                     # Partition
#SBATCH -n 1                             # Number of tasks (usually 1 for array jobs)
#SBATCH -c 2                             # CPUs per task
#SBATCH -t 4:00:00                            # Time (minutes)
#SBATCH -o /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/trim_testing/%x-%A_%a.log
#SBATCH -e /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/trim_testing/%x-%A_%a.log


INDEX=$SLURM_ARRAY_TASK_ID
# CD to directory
cd /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_Tag/common_files/data/external/spike_in

# Inform the app how many cores we requested for our job. The app can use this many cores.
# The special $NSLOTS keyword is automatically set to the number used on the -pe line above.
export OMP_NUM_THREADS=$NSLOTS

# activate all neeeded modules and packages
# source activate personal_software
module load functional_genomics/qc/fastqc/0.12.1
module load functional_genomics/tools/kentUtils/302
module load functional_genomics/tools/picard/3.0.0
activate_project /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/

bowtie2-inspect -n GRCh38_no_alt