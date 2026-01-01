import os
import logging
import subprocess
import glob


def run_sicer_peak_calling(Configuration):
    """
    Run SICER for CUT&RUN histone marks.
    Automatically activates the conda environment that contains SICER.
    """

    sample = Configuration.file_to_process

    # --- Skip SICER for IgG samples themselves ---
    if "IgG" in sample:
        logging.info(f"Skipping SICER for IgG sample: {sample}")
        return

    SICER_dir = Configuration.SICER_dir
    sample_SICER_dir = os.path.join(SICER_dir, sample)
    os.makedirs(sample_SICER_dir, exist_ok=True)

    alignment_dir = os.path.join(Configuration.Alignment_dir, sample)
    sicer_env = "sicer2-env"

    # -------------------------------------------------------------------------
    # Function to find matching IgG control BED
    # -------------------------------------------------------------------------
    def find_igg_control_bam(sample, alignment_root, tstatus):
        parts = sample.split("_")
        base_prefix = "_".join(parts[0:3])
        pattern = os.path.join(alignment_root, f"{base_prefix}_IgG")

        matches = glob.glob(pattern)
        if not matches:
            logging.warning(f"No IgG directory found for {sample}")
            return None

        igg_dir = matches[0]
        ctrl_bam_same = os.path.join(
            igg_dir, f"{os.path.basename(igg_dir)}.{tstatus}.test.bedgraph"
        )
        if os.path.exists(ctrl_bam_same):
            return ctrl_bam_same

        other_status = "rmDup" if tstatus == "dupMarked" else "dupMarked"
        ctrl_bam_other = os.path.join(
            igg_dir, f"{os.path.basename(igg_dir)}.{other_status}.test.bedgraph"
        )
        if os.path.exists(ctrl_bam_other):
            logging.warning(
                f"No IgG BED for {tstatus}, using {other_status} instead: {ctrl_bam_other}"
            )
            return ctrl_bam_other

        logging.warning(f"No IgG BED found for {sample} ({tstatus} or {other_status})")
        return None

    # -------------------------------------------------------------------------
    # Histone mark–specific SICER parameters
    # -------------------------------------------------------------------------
    if "H3K4me3" in sample or "H3K27ac" in sample or "H3K4me1" in sample:
        params = {"window_size": 200, "gap_size": 600, "evalue": 100}
    elif "H3K27me3" in sample:
        params = {"window_size": 3000, "gap_size": 30000, "evalue": 1000}
    elif "H3K36me3" in sample or "H3K9me3" in sample:
        params = {"window_size": 5000, "gap_size": 50000, "evalue": 5000}
    else:
        params = {"window_size": 200, "gap_size": 600, "evalue": 100}

    statuses = ["rmDup", "dupMarked"]
    results = {}

    for tstatus in statuses:
        treatment_bam = os.path.join(
            alignment_dir, f"{sample}.{tstatus}.test.bedgraph"
        )
        if not os.path.exists(treatment_bam):
            logging.warning(f"Treatment BED not found for {sample} ({tstatus}) — skipping.")
            continue

        control_bam = find_igg_control_bam(sample, Configuration.Alignment_dir, tstatus)

        # -------------------------------------------------------------------------
        # Convert bedGraph (4-column) → BED (6-column) with integer scores
        # -------------------------------------------------------------------------
        def convert_bedgraph_to_bed(input_bedgraph, output_bed, multiplier=1_000_000):
            """
            Convert a 4-column bedGraph (chrom, start, end, score)
            into a 6-column BED file (chrom, start, end, name, score, strand),
            scaling floating-point scores into integers.
            """
            with open(input_bedgraph, "r") as infile, open(output_bed, "w") as outfile:
                for line in infile:
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        chrom, start, end, score = parts[0], parts[1], parts[2], parts[3]
                        try:
                            score_int = int(round(float(score) * multiplier))
                        except ValueError:
                            score_int = 0
                        name = "."
                        strand = "."
                        outfile.write(
                            f"{chrom}\t{start}\t{end}\t{name}\t{score_int}\t{strand}\n"
                        )
            logging.info(
                f"Converted {input_bedgraph} to BED format with integer scores (×{multiplier}): {output_bed}"
            )

        # Convert treatment + control files
        treatment_bed = treatment_bam.replace(".bedgraph", ".bed")
        convert_bedgraph_to_bed(treatment_bam, treatment_bed)
        treatment_bam = treatment_bed

        if control_bam:
            control_bed = control_bam.replace(".bedgraph", ".bed")
            convert_bedgraph_to_bed(control_bam, control_bed)
            control_bam = control_bed

        # -------------------------------------------------------------------------
        # Build SICER command
        # -------------------------------------------------------------------------
        sicer_cmd = [
            "sicer",
            "-t", treatment_bam,
            "-s", "hg38",
            "-w", str(params["window_size"]),
            "-g", str(params["gap_size"]),
            "-e", str(params["evalue"]),
            "--significant_reads",
            "-o", sample_SICER_dir,
        ]

        if control_bam:
            sicer_cmd.insert(3, "-c")
            sicer_cmd.insert(4, control_bam)
            output_prefix = os.path.join(
                sample_SICER_dir, f"{sample}_treat-{tstatus}_ctrl"
            )
            logging.info(
                f"Running SICER for {sample} ({tstatus}) **with IgG control**: {control_bam}"
            )
        else:
            output_prefix = os.path.join(
                sample_SICER_dir, f"{sample}_treat-{tstatus}_noCtrl"
            )
            logging.warning(
                f"No IgG control found for {sample} ({tstatus}) — running SICER **without control**"
            )

        # -------------------------------------------------------------------------
        # Run SICER inside conda environment
        # -------------------------------------------------------------------------
        shell_cmd = (
            f"bash -c '"
            f"source $(conda info --base)/etc/profile.d/conda.sh && "
            f"conda activate {sicer_env} && "
            f"{' '.join(sicer_cmd)}'"
        )

        logging.info(f"Executing in env [{sicer_env}]: {shell_cmd}")

        try:
            subprocess.run(shell_cmd, shell=True, check=True)
            logging.info(f"SICER completed successfully for {sample} ({tstatus})")

            results[tstatus] = os.path.join(
                sample_SICER_dir,
                f"{os.path.basename(output_prefix)}-W{params['window_size']}-G{params['gap_size']}.scoreisland",
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"SICER failed for {sample} ({tstatus}) with error: {e}")
        except Exception as e:
            logging.error(f"Unexpected error running SICER for {sample} ({tstatus}): {e}")

    return results

