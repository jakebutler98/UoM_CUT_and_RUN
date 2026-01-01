import os
import glob
import logging
import subprocess

def FRiP(Configuration):
    """
    Calculate FRiP (Fraction of Reads in Peaks) for both rmDup and dupMarked BAMs,
    using both MACS2 and SEACR peak files if available.
    Skips gracefully if any file is missing.

    Output: <sample>_FRiP_summary.tsv
    """

    sample = Configuration.file_to_process
    alignment_dir = os.path.join(Configuration.Alignment_dir, sample)
    macs2_dir = os.path.join(Configuration.MACS2_dir, sample)
    seacr_dir = os.path.join(Configuration.SEACR_dir, sample)
    output_dir = Configuration.peak_qc_dir
    os.makedirs(output_dir, exist_ok=True)

    # Skip IgG samples
    if "IgG" in sample:
        logging.info(f"Skipping FRiP calculation for IgG sample: {sample}")
        return

    # Helper: count mapped reads
    def count_total_reads(bam):
        try:
            cmd = ["samtools", "view", "-c", "-F", "260", bam]
            return int(subprocess.check_output(cmd).strip())
        except subprocess.CalledProcessError:
            return 0

    # Helper: count reads overlapping peaks
    def count_reads_in_peaks(bam, peaks):
        try:
            intersect_cmd = ["bedtools", "intersect", "-u", "-a", bam, "-b", peaks]
            count_cmd = ["samtools", "view", "-c", "-"]
            p1 = subprocess.Popen(intersect_cmd, stdout=subprocess.PIPE)
            result = subprocess.check_output(count_cmd, stdin=p1.stdout)
            p1.stdout.close()
            return int(result.strip())
        except subprocess.CalledProcessError:
            return 0

    # Helper: find peak files for a peak caller
    def get_peak_files(peak_dir, caller):
        if not os.path.exists(peak_dir):
            return []
        if caller == "MACS2":
            patterns = ["*.narrowPeak", "*.broadPeak"]
        elif caller == "SEACR":
            # SEACR outputs .bed files, possibly with "relaxed" or "stringent" suffixes
            patterns = ["*.bed", "*.relaxed.bed", "*.stringent.bed"]
        else:
            patterns = ["*"]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(peak_dir, pattern)))
        return sorted(files)

    # Gather peak files for both peak callers
    macs2_peaks = get_peak_files(macs2_dir, "MACS2")
    seacr_peaks = get_peak_files(seacr_dir, "SEACR")

    if not macs2_peaks and not seacr_peaks:
        logging.warning(f"No MACS2 or SEACR peaks found for {sample}, skipping FRiP.")
        return None

    peak_sources = []
    if macs2_peaks:
        peak_sources.append(("MACS2", macs2_peaks[0]))
    if seacr_peaks:
        peak_sources.append(("SEACR", seacr_peaks[0]))

    results = []
    statuses = ["rmDup", "dupMarked"]

    for status in statuses:
        bam_path = os.path.join(alignment_dir, f"{sample}.{status}.mapped.sorted.bam")
        if not os.path.exists(bam_path):
            logging.warning(f"BAM not found for {sample} [{status}], skipping.")
            continue

        total_reads = count_total_reads(bam_path)
        if total_reads == 0:
            logging.warning(f"No reads counted for {sample} [{status}], skipping FRiP.")
            continue

        for peak_source, peak_file in peak_sources:
            if not os.path.exists(peak_file):
                logging.warning(f"{peak_source} peak file missing for {sample}, skipping.")
                continue

            reads_in_peaks = count_reads_in_peaks(bam_path, peak_file)
            frip_score = reads_in_peaks / total_reads if total_reads > 0 else 0

            results.append((sample, status, peak_source, total_reads, reads_in_peaks, frip_score))
            logging.info(
                f"{sample} [{status}] {peak_source} FRiP = {frip_score:.4f} "
                f"({reads_in_peaks}/{total_reads})"
            )

    # Write summary
    if not results:
        logging.warning(f"No FRiP results generated for {sample}.")
        return None

    out_file = os.path.join(output_dir, f"{sample}_FRiP_summary.tsv")
    with open(out_file, "w") as f:
        f.write("Sample\tStatus\tPeakCaller\tTotalReads\tReadsInPeaks\tFRiP\n")
        for row in results:
            f.write("\t".join(map(str, row)) + "\n")

    logging.info(f"FRiP summary written to: {out_file}")
    return results