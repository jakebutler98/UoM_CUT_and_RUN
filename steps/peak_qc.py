########################################
# script for peak qc
########################################

import os
import logging
import subprocess


def get_total_mapped_reads(flagstat_path):
    with open(flagstat_path, 'r') as f:
        for line in f:
            if "mapped" in line:
                return int(line.split()[0])
    raise ValueError("Mapped read count not found in flagstat file.")


def get_reads_in_peaks(bam_file, peak_file):
    intersect = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-a", bam_file, "-b", peak_file],
        stdout=subprocess.PIPE
    )
    output = subprocess.check_output(["samtools", "view", "-c", "-"], stdin=intersect.stdout)
    return int(output.decode().strip())


def FRiP(Configuration):
    sample = Configuration.file_to_process
    sample_dir = os.path.join(Configuration.Alignment_dir, sample)
    peak_qc_dir = os.path.join(Configuration.peak_qc_dir, sample)
    os.makedirs(peak_qc_dir, exist_ok=True)

    stats_file = os.path.join(sample_dir, "stats", f"{sample}.filtered.flagstat.txt")
    bam_file = os.path.join(sample_dir, f"{sample}.dupMarked.bam")

    total_mapped_reads = get_total_mapped_reads(stats_file)
    logging.info(f"[{sample}] Total mapped reads: {total_mapped_reads}")

    # Identify MACS2 peak file (narrow or broad)
    macs2_dir = os.path.join(Configuration.MACS2_dir, sample)
    macs2_peak_files = [
        os.path.join(macs2_dir, f"{sample}_treat-dupMarked_ctrl_dupMarked_macs2_peaks.narrowPeak"),
        os.path.join(macs2_dir, f"{sample}_treat-dupMarked_ctrl_dupMarked_macs2_peaks.broadPeak")
    ]
    macs2_peak_file = next((f for f in macs2_peak_files if os.path.exists(f)), None)

    if not macs2_peak_file:
        raise FileNotFoundError(f"No MACS2 peak file (narrow or broad) found for {sample}.")

    macs2_reads = get_reads_in_peaks(bam_file, macs2_peak_file)
    macs2_frip = macs2_reads / total_mapped_reads
    logging.info(f"[{sample}] MACS2 peaks: {macs2_reads} reads, FRiP = {macs2_frip:.4f}")

    # SEACR peaks
    seacr_peak_files = [
        os.path.join(Configuration.SEACR_dir, sample, f"{sample}_treat-dupMarked_ctrl_dupMarked_SEACR_non_relaxed.relaxed.bed"),
        os.path.join(Configuration.SEACR_dir, sample, f"{sample}_treat-dupMarked_ctrl_dupMarked_SEACR_non_stringent.stringent.bed"),
        os.path.join(Configuration.SEACR_dir, sample, f"{sample}_treat-dupMarked_NOCTRL_SEACR_top_relaxed.relaxed.bed"),
        os.path.join(Configuration.SEACR_dir, sample, f"{sample}_treat-dupMarked_NOCTRL_SEACR_top_stringent.stringent.bed")
    ]
    seacr_file = next((f for f in seacr_peak_files if os.path.exists(f)), None)

    if not os.path.exists(seacr_file):
        raise FileNotFoundError(f"SEACR peak file not found for {sample}.")

    seacr_reads = get_reads_in_peaks(bam_file, seacr_file)
    seacr_frip = seacr_reads / total_mapped_reads
    logging.info(f"[{sample}] SEACR peaks: {seacr_reads} reads, FRiP = {seacr_frip:.4f}")

    # Write summary
    summary_path = os.path.join(peak_qc_dir, f"{sample}_FRiP_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Sample: {sample}\n")
        f.write(f"Total mapped reads: {total_mapped_reads}\n\n")

        f.write(f"MACS2 peak file: {os.path.basename(macs2_peak_file)}\n")
        f.write(f"Reads in peaks: {macs2_reads}\n")
        f.write(f"FRiP (MACS2): {macs2_frip:.4f}\n\n")

        f.write(f"SEACR peak file: {os.path.basename(seacr_file)}\n")
        f.write(f"Reads in peaks: {seacr_reads}\n")
        f.write(f"FRiP (SEACR): {seacr_frip:.4f}\n")

    logging.info(f"[{sample}] FRiP summary written to {summary_path}")
