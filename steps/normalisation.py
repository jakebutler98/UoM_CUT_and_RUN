import os
import logging
import subprocess

def spike_in_normalization(Configuration):
    """
    Normalizes BAM file coverage using the spike-in scaling factor and generates a normalized BAM file.
    """

    # Ensure spike-in alignment was performed
    if Configuration.spike_in != "TRUE":
        logging.info("[info] Spike-in normalization skipped (spike-in alignment was not performed).")
        return

    # Define file paths
    spikein_log = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.spikein.bowtie2.log")
    input_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dedup.bam")
    output_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.normalized.bam")
    
    if not os.path.exists(spikein_log):
        logging.error(f"Spike-in alignment log file {spikein_log} not found!")
        raise Exception("Spike-in log file not found")

    if not os.path.exists(input_bam):
        logging.error(f"Sample BAM file {input_bam} not found!")
        raise Exception("Sample BAM file not found")

    # Extract total reads from spike-in log
    with open(spikein_log, "r") as log_file:
        log_content = log_file.read()
        total_reads = int([line.split()[0] for line in log_content.splitlines() if "reads; of these:" in line][0])
        align_ratio = float([line.split()[0] for line in log_content.splitlines() if "overall alignment" in line][0].replace("%", ""))
        spikein_reads = int(total_reads * align_ratio / 100)

    # Compute spike-in scaling factor
    if spikein_reads == 0:
        logging.warning("[warning] No spike-in reads detected. Normalization factor set to 1 (no scaling).")
        scaling_factor = 1.0
    else:
        scaling_factor = total_reads / spikein_reads
        logging.info(f"[info] Spike-in normalization factor for {Configuration.file_to_process}: {scaling_factor:.4f}")

    # Save scaling factor for reference
    scaling_factor_file = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.scaling_factor.txt")
    with open(scaling_factor_file, "w") as outfile:
        outfile.write(f"{scaling_factor:.4f}\n")
    
    logging.info(f"Spike-in normalization factor saved to {scaling_factor_file}")

    # Convert BAM to BedGraph using BEDTools
    temp_bedgraph = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.bedgraph")
    subprocess.run(["bedtools", "genomecov", "-ibam", input_bam, "-bg", "-scale", str(scaling_factor), "-g", Configuration.genome_chrom_sizes, "-o", temp_bedgraph], check=True)
    logging.info("Converted BAM to scaled BedGraph using BEDTools.")

    # Convert scaled BedGraph back to BAM
    temp_bed = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.bed")
    subprocess.run(["bedtools", "bedtobam", "-i", temp_bedgraph, "-g", Configuration.genome_chrom_sizes, "-o", temp_bed], check=True)
    logging.info("Converted scaled BedGraph to BED format.")

    subprocess.run(["samtools", "view", "-bS", temp_bed, "-o", output_bam], check=True)
    logging.info(f"Converted BED to normalized BAM: {output_bam}")

    # Index the normalized BAM file
    subprocess.run(["samtools", "index", output_bam], check=True)
    logging.info(f"Indexed normalized BAM file: {output_bam}.bai")

    # Cleanup temporary files
    os.remove(temp_bedgraph)
    os.remove(temp_bed)

    logging.info("Spike-in normalization completed successfully.")


