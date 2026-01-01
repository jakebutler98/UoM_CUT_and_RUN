#!/bin/bash --login
#SBATCH -J full_test_norm_alt                 # Job name
#SBATCH -p multicore                     # Partition
#SBATCH -n 1                             # Number of tasks (usually 1 for array jobs)
#SBATCH -c 8                             # CPUs per task
#SBATCH -t 6:00:00                            # Time (minutes)
#SBATCH -a 1-144                          # Array range
#SBATCH -o /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/full_test_23_10/%x-%A_%a.log
#SBATCH -e /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/full_test_23_10/%x-%A_%a.log


INDEX=$SLURM_ARRAY_TASK_ID
# CD to directory
cd /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing

# Inform the app how many cores we requested for our job. The app can use this many cores.
# The special $NSLOTS keyword is automatically set to the number used on the -pe line above.
export OMP_NUM_THREADS=$NSLOTS

# activate all neeeded modules and packages
# source activate personal_software
module load functional_genomics/qc/fastqc/0.12.1
module load functional_genomics/tools/kentUtils/302
module load functional_genomics/tools/picard/3.0.0
activate_project /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/
module load functional_genomics/fastq/fastp/0.23.4
module load functional_genomics/tools/bedtools/2.30.0
module load apps/sequencing/macs2/2.2.9.1
module load functional_genomics/mapping/bowtie2/2.5.1


# module load tools/java/1.8.0
SAMPLE=$(awk "NR==$INDEX" /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/scripts/samples.txt)



sleep $(($INDEX*20))
python ./scripts/main_CUT_RUN.py -i ${SAMPLE} \
  -s call_peaks \
  -s FRiP
