#!/bin/bash --login
#SBATCH -J insert_metrics
#SBATCH -t 1:00:00
#SBATCH -p multicore
#SBATCH -c 2
#SBATCH -n 1
#SBATCH --array=1-${N}    # we will substitute N dynamically
#SBATCH -o logs/%x_%A_%a.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ---- Inputs passed when submitting ----
HISTONE=$1        # e.g. H3K27ac
LISTFILE=$2       # e.g. H3K27ac_bams.txt

module load functional_genomics/tools/picard/3.0.0
activate_project /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/transcription_factors/

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LISTFILE")
sample=$(basename "$bam" .dupMarked.bam)

metrics_dir="/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/qc/insert_metrics/${HISTONE}"
mkdir -p "$metrics_dir"

picard CollectInsertSizeMetrics \
    -I "$bam" \
    -O "$metrics_dir/${sample}_insert_metrics.txt" \
    -M 0.5 \
    --Histogram_FILE "$metrics_dir/${sample}_insert_histogram.pdf"