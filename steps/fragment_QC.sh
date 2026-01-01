#!/bin/bash --login
#SBATCH -J multiQC_CUTRUN
#SBATCH -t 8:00:00
#SBATCH -p multicore
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -o /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/%x.o%j
#SBATCH -e /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_pipeline/logs/%x.o%j

# === ENVIRONMENT ===
activate_project /mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# === DIRECTORIES ===
bam_base="/mnt/jw01-aruk-home01/projects/psa_functional_genomics/Jake_Butler/OA_temp_storage/OA_CUTRUN_output/Alignments"
out_base="/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/multiBamSummary"
mkdir -p "$out_base"

echo "Running deepTools QC across 144 samples"
echo "BAM base: $bam_base"
echo "Output base: $out_base"
echo "Threads: $OMP_NUM_THREADS"

# === HELPER FUNCTION ===
run_qc_set() {
    local label="$1"
    shift
    local bam_files=("$@")
    local out_dir="${out_base}/${label}"
    mkdir -p "$out_dir"

    if [ "${#bam_files[@]}" -lt 2 ]; then
        echo "Skipping $label — fewer than 2 BAM files"
        return
    fi

    # Extract labels
    labels=()
    for bam in "${bam_files[@]}"; do
        labels+=("$(basename "${bam}" .dupMarked.bam)")
    done

    echo "=== Processing ${label} ==="
    echo "  Found ${#bam_files[@]} BAMs"

    # ensure all bams are indexed
    for bam in "${bam_files[@]}"; do
        if [ ! -f "${bam}.bai" ]; then
            echo "  Indexing BAM: ${bam}"
            samtools index "$bam"
        fi
    done

    # MultiBamSummary
    npz="$out_dir/${label}_multiBamSummary.npz"
    multiBamSummary bins -bs 500 \
        --bamfiles "${bam_files[@]}" \
        --labels "${labels[@]}" \
        --numberOfProcessors $OMP_NUM_THREADS \
        --outFileName "$npz" \
        -p $OMP_NUM_THREADS

    # PCA
    plotPCA \
        --corData "$npz" \
        --labels "${labels[@]}" \
        --plotFile "$out_dir/${label}_PCA.pdf" \
        --outFileNameData "$out_dir/${label}_PCA.tab"

    # Correlation
    plotCorrelation \
        -in "$npz" \
        --corMethod spearman --skipZeros \
        --labels "${labels[@]}" \
        --whatToPlot heatmap \
        --plotTitle "Spearman correlation - ${label}" \
        --plotFile "$out_dir/${label}_correlation_heatmap.pdf" \
        --outFileCorMatrix "$out_dir/${label}_correlation_matrix.tab"

    # Fingerprint (optional; can be large)
    plotFingerprint \
        --bamfiles "${bam_files[@]}" \
        --labels "${labels[@]}" \
        --plotFile "$out_dir/${label}_fingerprint.pdf" \
        --outRawCounts "$out_dir/${label}_fingerprint_counts.txt" \
        --outQualityMetrics "$out_dir/${label}_fingerprint_metrics.txt" \
        --numberOfProcessors $OMP_NUM_THREADS

    echo "=== Finished ${label} ==="
}

# === STEP 1: Detect histone marks ===
marks=$(ls "$bam_base" | grep -oE 'H3K[0-9A-Za-z]+(me[0-9])?|CTCF|IgG' | sort -u)
conditions=("H" "D")

# === STEP 2: Run comparisons ===

## (A) Global — all BAMs together
#mapfile -t all_bams < <(find "$bam_base" -type f -name "*.dupMarked.bam" | sort)
#run_qc_set "GLOBAL_all_marks" "${all_bams[@]}"

## (B) Per histone mark (all conditions together)
for mark in $marks; do
    echo "Processing mark: $mark"
    mapfile -t mark_bams < <(find "$bam_base" -type f -path "*${mark}*" -name "*.dupMarked.bam" | sort)
    run_qc_set "${mark}_ALL" "${mark_bams[@]}"

    ## (C) Per mark × condition
    for cond in "${conditions[@]}"; do
        mapfile -t cond_bams < <(find "$bam_base" -type f -path "*_C_${cond}_*" -path "*${mark}*" -name "*.dupMarked.bam" | sort)
        if [ "${#cond_bams[@]}" -ge 2 ]; then
            run_qc_set "${mark}_${cond}" "${cond_bams[@]}"
        fi
    done
done

echo "All comparisons completed."