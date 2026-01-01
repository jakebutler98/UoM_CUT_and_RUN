########################################
# script for peak calling


########################################

import os
import logging
import subprocess
import glob

def run_macs2_peak_calling(Configuration):
    """
    Runs MACS2 peak calling for both rmDup and dupMarked BAMs.
    Flexible IgG control matching (suffix, status, fallback).
    If IgG control is missing, runs MACS2 without control.
    Skips MACS2 for IgG samples.
    """

    def find_igg_control_bam(sample, alignment_dir, tstatus):
        parts = sample.split("_")
        base_prefix = "_".join(parts[0:3])  # e.g., OA004_C_H
        read_suffix = parts[-1] if parts[-1] in ("R1", "R2") else None

        # 1. Prefer IgG with same read suffix
        if read_suffix:
            pattern_same = os.path.join(alignment_dir, f"{base_prefix}_IgG_{read_suffix}")
            matches_same = glob.glob(pattern_same)
            if matches_same:
                ctrl_dir = matches_same[0]
                ctrl_file_same = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{tstatus}.mapped.sorted.bam")
                if os.path.exists(ctrl_file_same):
                    return ctrl_file_same, tstatus
                # fallback to opposite status if same not found
                other_status = "rmDup" if tstatus == "dupMarked" else "dupMarked"
                ctrl_file_other = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{other_status}.mapped.sorted.bam")
                if os.path.exists(ctrl_file_other):
                    return ctrl_file_other, other_status

        # 2. Fall back to IgG without read suffix
        pattern_no_read = os.path.join(alignment_dir, f"{base_prefix}_IgG")
        matches_no_read = glob.glob(pattern_no_read)
        if matches_no_read:
            ctrl_dir = matches_no_read[0]
            ctrl_file_same = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{tstatus}.mapped.sorted.bam")
            if os.path.exists(ctrl_file_same):
                return ctrl_file_same, tstatus
            # fallback to opposite status if same not found
            other_status = "rmDup" if tstatus == "dupMarked" else "dupMarked"
            ctrl_file_other = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{other_status}.mapped.sorted.bam")
            if os.path.exists(ctrl_file_other):
                return ctrl_file_other, other_status

        return None, None

    sample = Configuration.file_to_process

    # --- Skip MACS2 for IgG samples ---
    if "IgG" in sample:
        logging.info(f"Skipping MACS2 for IgG sample: {sample}")
        return

    MACS2_dir = Configuration.MACS2_dir
    sample_MACS2_dir = os.path.join(MACS2_dir, sample)
    os.makedirs(sample_MACS2_dir, exist_ok=True)
    alignment_dir = os.path.join(Configuration.Alignment_dir, sample)
    statuses = ["rmDup", "dupMarked"]

    def get_peak_type(sample):
        is_broad = any(mark in sample for mark in ["H3K27me3", "H3K36me3", "H3K9me3"])
        is_narrow = any(mark in sample for mark in ["H3K4me3", "H3K27ac", "H3K4me1"])
        is_tf = not (sample.upper().startswith("H3") or sample.upper().startswith("H4"))
        if is_tf:
            is_narrow = True
        return is_broad, is_narrow

    is_broad, is_narrow = get_peak_type(sample)

    results = {}
    for tstatus in statuses:
        treatment_bam = os.path.join(alignment_dir, f"{sample}.{tstatus}.mapped.sorted.bam")
        if not os.path.exists(treatment_bam):
            logging.warning(f"Treatment BAM not found for {sample} ({tstatus}) — skipping.")
            continue

        # Try to find an IgG control BAM
        control_bam, cstatus = find_igg_control_bam(sample, Configuration.Alignment_dir, tstatus)

        macs2_output_prefix = os.path.join(
            sample_MACS2_dir,
            f"{sample}_treat-{tstatus}" + (f"_ctrl_{cstatus}" if cstatus else "") + "_macs2"
        )

        macs2_cmd = [
            "macs2", "callpeak",
            "-t", treatment_bam,
            "--outdir", sample_MACS2_dir,
            "-n", macs2_output_prefix,
            "--format", "BAMPE",
            "--gsize", "hs",
            "--nomodel", "--keep-dup", "all"
        ]

        if control_bam:
            macs2_cmd.extend(["-c", control_bam])
            logging.info(f"Running MACS2: {os.path.basename(treatment_bam)} vs {os.path.basename(control_bam)}")
        else:
            logging.warning(f"No IgG control BAM exists for {sample} ({tstatus}) — running MACS2 without control.")
            logging.info(f"Running MACS2: {os.path.basename(treatment_bam)} [no control]")

        if is_broad:
            macs2_cmd.extend(["--broad", "--broad-cutoff", "0.1"])

        # --- Run MACS2 inside the conda environment ---
        macs2_run_cmd = ["conda", "run", "-n", "macs2-env", "bash", "-c", " ".join(macs2_cmd)]
        logging.debug(f"Executing: {' '.join(macs2_run_cmd)}")
        subprocess.run(macs2_run_cmd, check=True)

        logging.info(f"MACS2 peak calling complete for {tstatus}")

        # Output file
        results[tstatus] = (
            os.path.join(
                sample_MACS2_dir,
                f"{sample}_treat-{tstatus}" + (f"_ctrl_{cstatus}" if cstatus else "")
                + ("_macs2_peaks.broadPeak" if is_broad else "_macs2_peaks.narrowPeak")
            )
        )

    logging.info("MACS2 peak calling completed for both rmDup and dupMarked.")
    return results





### DEVELOPMENT PEAK CALLING FUNCTION ###

def run_peak_calling_v2(Configuration):
    """
    Runs SEACR peak calling with:
      - Threshold policy based on histone mark
      - Flexible IgG control matching:
            * Prefer IgG with same read suffix (_R1/_R2)
            * Fall back to IgG without suffix
            * Fall back to opposite status if same not found
      - If IgG missing:
            * If IgG exists but not processed, skip until ready
            * If IgG does not exist, run SEACR top % peaks mode
      - If sample is IgG, skip SEACR
    """

    def threshold_for_mark(mark):
        if mark in ["H3K27ac", "H3K4me1", "H3K4me3"]:
            return "relaxed"
        elif mark in ["H3K27me3", "H3K36me3", "H3K9me3"]:
            return "stringent"
        else:
            return "relaxed"

    def find_igg_control(sample, alignment_dir, tstatus):
        parts = sample.split("_")
        base_prefix = "_".join(parts[0:3])  # e.g., OA004_C_H
        read_suffix = parts[-1] if parts[-1] in ("R1", "R2") else None

        # 1. Prefer IgG with same read suffix
        if read_suffix:
            pattern_same = os.path.join(alignment_dir, f"{base_prefix}_IgG_{read_suffix}")
            matches_same = glob.glob(pattern_same)
            if matches_same:
                ctrl_dir = matches_same[0]
                ctrl_file_same = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{tstatus}.bowtie2.fragments.normalised.bedgraph")
                if os.path.exists(ctrl_file_same):
                    return ctrl_file_same, tstatus
                # fallback to opposite status if same not found
                other_status = "rmDup" if tstatus == "dupMarked" else "dupMarked"
                ctrl_file_other = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{other_status}.bowtie2.fragments.normalised.bedgraph")
                if os.path.exists(ctrl_file_other):
                    return ctrl_file_other, other_status

        # 2. Fall back to IgG without read suffix
        pattern_no_read = os.path.join(alignment_dir, f"{base_prefix}_IgG")
        matches_no_read = glob.glob(pattern_no_read)
        if matches_no_read:
            ctrl_dir = matches_no_read[0]
            ctrl_file_same = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{tstatus}.bowtie2.fragments.normalised.bedgraph")
            if os.path.exists(ctrl_file_same):
                return ctrl_file_same, tstatus
            # fallback to opposite status if same not found
            other_status = "rmDup" if tstatus == "dupMarked" else "dupMarked"
            ctrl_file_other = os.path.join(ctrl_dir, f"{os.path.basename(ctrl_dir)}.{other_status}.bowtie2.fragments.normalised.bedgraph")
            if os.path.exists(ctrl_file_other):
                return ctrl_file_other, other_status

        return None, None

    sample = Configuration.file_to_process

    # --- Skip SEACR for IgG samples ---
    if "IgG" in sample:
        logging.info(f"Skipping SEACR for IgG sample: {sample}")
        return

    SEACR_dir = Configuration.SEACR_dir
    sample_SEACR_dir = os.path.join(SEACR_dir, sample)
    os.makedirs(sample_SEACR_dir, exist_ok=True)

    # Parse histone mark
    try:
        mark = sample.split("_")[3]  # e.g., OA14_C_D_H3K27ac -> H3K27ac
    except IndexError:
        logging.warning(f"Could not parse mark from sample name: {sample}")
        mark = ""

    seacr_threshold = threshold_for_mark(mark)
    logging.info(f"SEACR threshold for {mark}: {seacr_threshold}")

    bedgraph_dir = os.path.join(Configuration.Alignment_dir, sample)
    statuses = ["dupMarked", "rmDup"]
    seacr_script_path = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/envs/default/bin/SEACR_1.3.sh"

    for tstatus in statuses:
        treatment_file = os.path.join(bedgraph_dir, f"{sample}.{tstatus}.bowtie2.fragments.normalised.bedgraph")
        if not os.path.exists(treatment_file):
            logging.warning(f"Treatment file not found for {sample} ({tstatus}) — skipping.")
            continue

        # Try to find an IgG control
        control_file, cstatus = find_igg_control(sample, Configuration.Alignment_dir, tstatus)

        if not control_file:
            # --- Check if IgG folder exists (unprocessed case) ---
            parts = sample.split("_")
            base_prefix = "_".join(parts[0:3])
            read_suffix = parts[-1] if parts[-1] in ("R1", "R2") else None

            raw_patterns = []
            if read_suffix:
                raw_patterns.append(os.path.join(Configuration.Alignment_dir, f"{base_prefix}_IgG_{read_suffix}"))
            raw_patterns.append(os.path.join(Configuration.Alignment_dir, f"{base_prefix}_IgG"))

            igg_exists_but_unprocessed = False
            for p in raw_patterns:
                matches = glob.glob(p)
                if matches:
                    igg_exists_but_unprocessed = True
                    logging.warning(
                        f"IgG control folder found for {sample} but no {tstatus}.bowtie2.fragments.normalised.bedgraph yet. "
                        "Likely still processing — skipping SEACR for now."
                    )
                    break

            if igg_exists_but_unprocessed:
                continue  # Wait until IgG is ready

            # --- No IgG exists → run top % peaks mode ---
            logging.warning(f"No IgG control exists for {sample} — running SEACR top peaks mode.")
            seacr_output_prefix = os.path.join(
                sample_SEACR_dir,
                f"{sample}_treat-{tstatus}_NOCTRL_SEACR_top_{seacr_threshold}"
            )
            seacr_cmd = [
                "conda", "run", "-n", "seacr_env",
                "bash", seacr_script_path,
                treatment_file,
                "0.01",  # top 1% peaks — adjust if needed
                "non",
                seacr_threshold,
                seacr_output_prefix
            ]
            logging.info(f"Running SEACR (top peaks mode): {os.path.basename(treatment_file)} [non, {seacr_threshold}]")
            subprocess.run(seacr_cmd, check=True)
            logging.info(f"SEACR top peaks mode complete. Output: {seacr_output_prefix}")
            continue

        # --- Standard SEACR run with IgG control ---
        seacr_output_prefix = os.path.join(
            sample_SEACR_dir,
            f"{sample}_treat-{tstatus}_ctrl_{cstatus}_SEACR_non_{seacr_threshold}"
        )
        seacr_cmd = [
            "conda", "run", "-n", "seacr_env",
            "bash", seacr_script_path,
            treatment_file,
            control_file,
            "non",
            seacr_threshold,
            seacr_output_prefix
        ]
        logging.info(f"Running SEACR: {os.path.basename(treatment_file)} vs {os.path.basename(control_file)} "
                     f"[non, {seacr_threshold}]")
        subprocess.run(seacr_cmd, check=True)
        logging.info(f"SEACR peak calling complete. Output: {seacr_output_prefix}")

