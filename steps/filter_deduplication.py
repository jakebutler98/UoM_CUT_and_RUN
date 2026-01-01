import os
import logging
import subprocess
import re

def filter_reads(Configuration):
    input_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}_sorted.bam")
    filtered_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.filtered.bam")

    if not os.path.exists(input_bam):
        logging.error(f"Input BAM file {input_bam} not found!")
        raise Exception("BAM file not found")

    logging.info(f"Starting read filtering for {Configuration.file_to_process}")

    step1_bam = filtered_bam.replace(".bam", ".step1.bam")
    step2_bam = filtered_bam.replace(".bam", ".step2.bam")
    step3_bam = filtered_bam.replace(".bam", ".step3.bam")
    step4_bam = filtered_bam.replace(".bam", ".step4.bam")

    # filter unmapped and non primary alignments - step 1
    subprocess.run(["samtools", "view", "-F", "260", "-f", "2", "-b", "-@", Configuration.num_cores, input_bam, "-o", step1_bam], check=True)
    
    # quality thresholding (hardcoded at 10) - step 2
    subprocess.run(["samtools", "view", "-q", str(Configuration.mapping_quality_threshold), "-b", "-@", Configuration.num_cores, step1_bam, "-o", step2_bam], check=True)
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step2_bam], check=True)

    # remove blacklisted - step 3
    if not os.path.exists(Configuration.blacklist_bed):
        logging.error(f"Blacklist BED file not found: {Configuration.blacklist_bed}")
        raise Exception("Blacklist BED file missing")

    with open(step3_bam, "wb") as outfile:
        subprocess.run(["bedtools", "intersect", "-v", "-a", step2_bam, "-b", Configuration.blacklist_bed], stdout=outfile, check=True)
    if not os.path.exists(step3_bam):
        logging.error(f"Step 3 BAM not created: {step3_bam}. bedtools may have failed.")
        raise Exception("Step 3 BAM file not generated")
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step3_bam], check=True)
    result = subprocess.run(["samtools", "idxstats", "-@", Configuration.num_cores, step3_bam], capture_output=True, text=True, check=True)
    
    # Filter out mitochondrial contigs - step 4
    contigs = [line.split('\t')[0] for line in result.stdout.strip().split('\n') if not line.startswith("chrM")]
    with open(step4_bam, "wb") as outfile:
        subprocess.run(["samtools", "view", "-b", "-@", Configuration.num_cores, step3_bam] + contigs, stdout=outfile, check=True)
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step4_bam], check=True)

    # Move final filtered file to its intended name
    os.rename(step4_bam, filtered_bam)
    os.rename(step4_bam + ".bai", filtered_bam + ".bai")

    # Cleanup intermediate files
    for bam in [step1_bam, step2_bam, step3_bam]:
        bai = bam + ".bai"
        try:
            os.remove(bam)
            if os.path.exists(bai):
                os.remove(bai)
        except Exception as e:
            logging.warning(f"Failed to remove intermediate file {bam}: {e}")

    logging.info(f"Read filtering completed. Filtered BAM saved as {filtered_bam}")
    
    
    # collect flagstats on filtered bam for future normalisation
    stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")
    filtered_flagstats = os.path.join(stats_dir, f"{Configuration.file_to_process}.filtered.flagstat.txt")
    
    os.makedirs(stats_dir, exist_ok=True)
    with open(filtered_flagstats, "w") as outfile:
        subprocess.run(["samtools", "flagstat", filtered_bam], stdout=outfile, check=True)
    logging.info(f"Flagstats for filtered BAM saved to {filtered_flagstats}")
    logging.info("Filtered BAM flagstats collected.")

    # delete the input BAM file to save space
    try:
        os.remove(input_bam)
        logging.info(f"Deleted input BAM file to save space: {input_bam}")
    except Exception as e:
        logging.warning(f"Failed to remove input BAM file {input_bam}: {e}")

    
    return filtered_bam


def spike_in_filter(Configuration):
    """Filter out reads that do NOT align uniquely."""
    alignment_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process)
    spikein_dir = os.path.join(alignment_dir, "spikein")
    input_bam = os.path.join(spikein_dir, f"{Configuration.file_to_process}_spikein_sorted.bam")
    filtered_bam = os.path.join(spikein_dir, f"{Configuration.file_to_process}_spikein_filtered.bam")
    
    if not os.path.exists(input_bam):
        logging.error(f"Input BAM file {input_bam} not found!")
        raise Exception("BAM file not found")
    
    logging.info(f"Starting spike-in filtering for {Configuration.file_to_process}")
    
    step1_bam = filtered_bam.replace(".bam", ".step1.bam")
    step2_bam = filtered_bam.replace(".bam", ".step2.bam")
    step3_bam = filtered_bam.replace(".bam", ".step3.bam")
    step4_bam = filtered_bam.replace(".bam", ".step4.bam")
    
    # filter unmapped and non primary alignments - step 1
    subprocess.run(["samtools", "view", "-q", str(Configuration.mapping_quality_threshold), "-F", "260", "-f", "2", "-b", "-@", Configuration.num_cores, input_bam, "-o", step1_bam], check=True)
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step1_bam], check=True)
    
    # quality thresholding (hardcoded at 10) - step 2
    subprocess.run(["samtools", "view", "-q", str(Configuration.mapping_quality_threshold), "-b", "-@", Configuration.num_cores, step1_bam, "-o", step2_bam], check=True)
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step2_bam], check=True)

    # remove blacklisted - step 3
    if not os.path.exists(Configuration.blacklist_bed):
        logging.error(f"Blacklist BED file not found: {Configuration.blacklist_bed}")
        raise Exception("Blacklist BED file missing")

    with open(step3_bam, "wb") as outfile:
        subprocess.run(["bedtools", "intersect", "-v", "-a", step2_bam, "-b", Configuration.blacklist_bed], stdout=outfile, check=True)
    if not os.path.exists(step3_bam):
        logging.error(f"Step 3 BAM not created: {step3_bam}. bedtools may have failed.")
        raise Exception("Step 3 BAM file not generated")
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step3_bam], check=True)
    result = subprocess.run(["samtools", "idxstats", "-@", Configuration.num_cores, step3_bam], capture_output=True, text=True, check=True)
    
    # Filter out mitochondrial contigs - step 4
    contigs = [line.split('\t')[0] for line in result.stdout.strip().split('\n') if not line.startswith("chrM")]
    with open(step4_bam, "wb") as outfile:
        subprocess.run(["samtools", "view", "-b", "-@", Configuration.num_cores, step3_bam] + contigs, stdout=outfile, check=True)
    subprocess.run(["samtools", "index", "-@", Configuration.num_cores, step4_bam], check=True)

    # Move final filtered file to its intended name
    os.rename(step4_bam, filtered_bam)
    os.rename(step4_bam + ".bai", filtered_bam + ".bai")

    # Cleanup intermediate files
    for bam in [step1_bam, step2_bam, step3_bam]:
        bai = bam + ".bai"
        try:
            os.remove(bam)
            if os.path.exists(bai):
                os.remove(bai)
        except Exception as e:
            logging.warning(f"Failed to remove intermediate file {bam}: {e}")
    
    logging.info(f"Spike-in filtering completed. Filtered BAM saved as {filtered_bam}")
    
    # calculate flagstats for filtered spike-in
    stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")
    filtered_flagstats = os.path.join(stats_dir, f"{Configuration.file_to_process}.spikein.filtered.flagstat.txt")
    os.makedirs(stats_dir, exist_ok=True)
    
    with open(filtered_flagstats, "w") as outfile:
        subprocess.run(["samtools", "flagstat", filtered_bam], stdout=outfile, check=True)
    logging.info(f"Flagstats for filtered spike-in BAM saved to {filtered_flagstats}")
    
    return filtered_bam



def mark_and_remove_duplicates(Configuration):
    """
    Marks and removes duplicates from the filtered BAM file using Picard's MarkDuplicates.
    """

    # Paths
    filtered_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.filtered.bam")
    dupMarked_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dupMarked.bam")
    metrics_file_dupMark = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dupMark.txt")
    rmDup_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.rmDup.bam")
    metrics_file_rmDup = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.rmDup.txt")


    if not os.path.exists(filtered_bam):
        logging.error(f"Filtered BAM file {filtered_bam} not found for deduplication!")
        raise Exception("Filtered BAM not found")

    logging.info(f"Marking duplicates for {Configuration.file_to_process}")

    # Run Picard MarkDuplicates - for dupMarking
    picard_command = [
        "picard", "MarkDuplicates",
        f"I={filtered_bam}",
        f"O={dupMarked_bam}",
        f"M={metrics_file_dupMark}",
        f"REMOVE_DUPLICATES=false",
    ]

    subprocess.run(picard_command)

    # Run Picard MarkDuplicates - for deduplication
    picard_command = [
        "picard", "MarkDuplicates",
        f"I={filtered_bam}",
        f"O={rmDup_bam}",
        f"M={metrics_file_rmDup}",
        "REMOVE_DUPLICATES=true"
    ]
    subprocess.run(picard_command)

    # Index final BAM
    subprocess.run(["samtools", "index", dupMarked_bam])
    subprocess.run(["samtools", "index", rmDup_bam])

    logging.info(f"Duplicate marking complete. Final BAM: {dupMarked_bam}")
    logging.info(f"Duplicate removal complete. Final BAM: {rmDup_bam}")


    ### FRAGMENT SIZE DISTRIBUTION

    fragmentLength_file_dupMark = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dupMark.fragmentLen.txt")
    fragmentLength_file_rmDup = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.rmDup.fragmentLen.txt")

    # Command for dupMarked BAM
    fragment_length_command_dupMark = (
        f"samtools view -F 0x04 {dupMarked_bam} | "
        "awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | "
        "sort | uniq -c | "
        f"awk -v OFS='\t' '{{print $2, $1/2}}' > {fragmentLength_file_dupMark}"
    )
    subprocess.run(fragment_length_command_dupMark, shell=True, check=True)

    # Command for rmDup BAM
    fragment_length_command_rmDup = (
        f"samtools view -F 0x04 {rmDup_bam} | "
        "awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | "
        "sort | uniq -c | "
        f"awk -v OFS='\t' '{{print $2, $1/2}}' > {fragmentLength_file_rmDup}"
    )
    subprocess.run(fragment_length_command_rmDup, shell=True, check=True)

    logging.info(f"Fragment length distribution saved to {fragmentLength_file_dupMark} and {fragmentLength_file_rmDup}")

    # Can remove the filtered BAM file to save space
    #try:
        #os.remove(filtered_bam)
        #logging.info(f"Deleted filtered BAM file to save space: {filtered_bam}")
    #except Exception as e:
        #logging.warning(f"Failed to remove filtered BAM file {filtered_bam}: {e}")
    

    return dupMarked_bam, rmDup_bam, metrics_file_dupMark, metrics_file_rmDup, fragmentLength_file_dupMark, fragmentLength_file_rmDup


def calculate_normalisation(Configuration):
    """
    Calculates spike-in normalization factor based on filtered BAM files.
    """

    logging.info("Starting spike-in normalization factor calculation...")

    stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")
    spikein_flagstat = os.path.join(stats_dir, f"{Configuration.file_to_process}.spikein.filtered.flagstat.txt")
    filtered_flagstat = os.path.join(stats_dir, f"{Configuration.file_to_process}.filtered.flagstat.txt")

    if not os.path.exists(spikein_flagstat) or not os.path.exists(filtered_flagstat):
        logging.error("One or both flagstat files are missing.")
        raise Exception("Missing flagstat files for normalization.")

    def get_mapped_reads(flagstat_path):
        with open(flagstat_path, "r") as f:
            for line in f:
                if " mapped (" in line and not line.startswith("0 "):
                    # Extract the first number (mapped reads)
                    match = re.match(r"(\d+)", line)
                    if match:
                        return int(match.group(1))
        raise Exception(f"Could not parse mapped reads from {flagstat_path}")

    # Extract mapped reads
    spikein_reads = get_mapped_reads(spikein_flagstat)
    experimental_reads = get_mapped_reads(filtered_flagstat)

    logging.info(f"Spike-in mapped reads: {spikein_reads}")
    logging.info(f"Experimental mapped reads: {experimental_reads}")

    if spikein_reads == 0:
        logging.error("Spike-in mapped reads is 0. Cannot compute normalization factor.")
        raise Exception("Invalid spike-in read count.")

    # Calculate normalization factor
    C = 10000  # Example constant
    norm_factor = round(C / spikein_reads, 6)

    # Save to a text file for downstream tools
    norm_file = os.path.join(stats_dir, f"{Configuration.file_to_process}.normalisation.txt")
    with open(norm_file, "w") as out:
        out.write(f"Spike-in mapped reads: {spikein_reads}\n")
        out.write(f"Experimental mapped reads: {experimental_reads}\n")
        out.write(f"Normalization factor: {norm_factor:.6f}\n")

    logging.info(f"Normalization factor calculated: {norm_factor}")
    logging.info(f"Saved to {norm_file}")

    return norm_factor


