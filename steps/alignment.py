import os
import logging
import subprocess


def get_sequence_length(Configuration):
    """
    Extracts the sequence length from FastQC output for the trimmed reads.
    """
    # Define the QC directory for the trimmed files
    qc_dir = os.path.join(Configuration.Trimmed_QC_dir, Configuration.file_to_process)

    # Check if the QC directory exists
    if not os.path.exists(qc_dir):
        raise Exception(f"QC directory {qc_dir} not found. Ensure QC has been run.")

    # Identify the FastQC zip file for R1 (assuming both R1 and R2 exist, but we take one)
    fastqc_zip = None
    for file in os.listdir(qc_dir):
        if file.endswith("1_fastqc.zip"):  # Checking only R1
            fastqc_zip = os.path.join(qc_dir, file)
            break

    if not fastqc_zip:
        raise Exception(f"No FastQC zip file found in {qc_dir}")

    # Unzip the FastQC output
    subprocess.run(["unzip", "-o", fastqc_zip, "-d", qc_dir], check=True)

    # Locate the extracted FastQC folder (same name as zip minus .zip)
    fastqc_folder = fastqc_zip.replace(".zip", "")

    # Path to FastQC data file
    fastqc_data_file = os.path.join(fastqc_folder, "fastqc_data.txt")

    # Ensure the data file exists
    if not os.path.exists(fastqc_data_file):
        raise Exception(f"FastQC data file not found: {fastqc_data_file}")

    # Extract sequence length
    with open(fastqc_data_file, "r") as file:
        for line in file:
            if line.startswith("Sequence length"):
                seq_length = line.split()[-1]
                # Handle cases where seq_length is a range (e.g. 0-76)
                if '-' in seq_length:
                    seq_length = int(seq_length.split('-')[-1])
                else:
                    seq_length = int(seq_length) 
                return seq_length  # Return the extracted length

    raise Exception("Sequence length not found in FastQC data")



def run_bowtie2_alignment(Configuration):
    """
    Runs Bowtie2 to align trimmed FASTQ files to a reference genome.
    """

    # Get the sequence length using the separate function
    seq_length = get_sequence_length(Configuration)
    logging.info(f"Extracted sequence length from FastQC: {seq_length} bp")

    # Set frag_120 based on sequence length
    Configuration.frag_120 = "TRUE" if seq_length <= 120 else "FALSE"

     # Define directories
    trimmed_dir = os.path.join(Configuration.Trimmed_dir, Configuration.file_to_process)
    alignment_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process)
    spikein_dir = os.path.join(alignment_dir, "spikein")  # Directory for spike-in alignment
    os.makedirs(alignment_dir, exist_ok=True)
    os.makedirs(spikein_dir, exist_ok=True)


    # Paths to trimmed paired-end files
    trim_R1 = os.path.join(trimmed_dir, f"{Configuration.file_to_process}_R1.fastq.gz")
    trim_R2 = os.path.join(trimmed_dir, f"{Configuration.file_to_process}_R2.fastq.gz")

    # Check if trimmed files exist
    if not os.path.exists(trim_R1) or not os.path.exists(trim_R2):
        logging.error(f"Trimmed files not found for {Configuration.file_to_process}. Please run trimming first.")
        raise Exception("Trimmed files not found")

    # Reference genome index path
    ref_genome_index = os.path.join(Configuration.Bowtie2_idx_dir)

    # Select Bowtie2 parameters based on frag_120
    if Configuration.frag_120 == "TRUE":
        logging.info("[info] Using dovetail mode for paired-end reads (frag_120 is TRUE).")
        bowtie2_command = [
            "bowtie2", "-p", str(Configuration.num_cores), "--dovetail", "--phred33", 
            "-x", ref_genome_index, "-1", trim_R1, "-2", trim_R2
        ]
    else:
        logging.info("[info] Using very-sensitive mode (frag_120 is FALSE).")
        bowtie2_command = [
            "bowtie2", "-p", str(Configuration.num_cores), "--end-to-end", "--very-sensitive", "--no-mixed", "--no-discordant", "--phred33", 
            "-I", "10", "-X", "700", "-x", ref_genome_index, "-1", trim_R1, "-2", trim_R2
        ]


    # Define output paths
    sam_output = os.path.join(alignment_dir, f"{Configuration.file_to_process}.sam")
    sorted_bam = os.path.join(alignment_dir, f"{Configuration.file_to_process}_sorted.bam")
    log_file_path = os.path.join(alignment_dir, f"{Configuration.file_to_process}.bowtie2.log")

    # Run Bowtie2 alignment and convert SAM to BAM, then sort
    with open(log_file_path, "w") as log_file:
        aligner = subprocess.Popen(bowtie2_command, stdout=subprocess.PIPE, stderr=log_file)
        converter = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=aligner.stdout, stdout=subprocess.PIPE, stderr=log_file)
        sorter = subprocess.Popen(["samtools", "sort", "-o", sorted_bam, "-"], stdin=converter.stdout, stderr=log_file)
        sorter.wait()


    # delete the SAM file to save space
    if os.path.exists(sam_output):
        os.remove(sam_output)
        logging.info(f"Deleted SAM file to save space: {sam_output}")

    # Index the sorted BAM file
    subprocess.run(["samtools", "index", sorted_bam])

    logging.info(f"Bowtie2 alignment complete for {Configuration.file_to_process}. Sorted BAM file saved to {sorted_bam}.")
    
    # Align to the spike-in reference genome
    if Configuration.spike_in:
        spikein_bam = os.path.join(spikein_dir, f"{Configuration.file_to_process}_spikein_sorted.bam")
        spikein_log = os.path.join(spikein_dir, f"{Configuration.file_to_process}.spikein.bowtie2.log")
        spikein_stats_file = os.path.join(spikein_dir, f"{Configuration.file_to_process}_spikein_stats.txt")

        logging.info(f"Running Bowtie2 alignment for {Configuration.file_to_process} against spike-in genome.")

        with open(spikein_log, "w") as log_file:
            aligner = subprocess.Popen(["bowtie2", "-p", str(Configuration.num_cores),"--end-to-end", "--very-sensitive", "--no-dovetail",  
            "--no-overlap", "--no-mixed", "--phred33", "--no-discordant", "-I", "10", "-X", "700",
            "-x", Configuration.spikein_genome, "-1", trim_R1, "-2", trim_R2],
                                        stdout=subprocess.PIPE, stderr=log_file)
            converter = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=aligner.stdout, stdout=subprocess.PIPE, stderr=log_file)
            sorter = subprocess.Popen(["samtools", "sort", "-o", spikein_bam, "-"], stdin=converter.stdout, stderr=log_file)
            sorter.wait()

        # Index the sorted BAM file for spike-in
        subprocess.run(["samtools", "index", spikein_bam])

    # Can now remove the trimmed FASTQ files to save space
    if os.path.exists(trim_R1):
        os.remove(trim_R1)
        logging.info(f"Deleted trimmed R1 FASTQ file to save space: {trim_R1}")
    if os.path.exists(trim_R2):
        os.remove(trim_R2)
        logging.info(f"Deleted trimmed R2 FASTQ file to save space: {trim_R2}")

    else:
        logging.info("[info] FASTQ files won't be aligned to the spike-in genome.")




def generate_alignment_stats(Configuration):
    """
    Generates alignment statistics using Samtools for the aligned BAM files.
    Saves statistics to log files in the stats directory.
    """
    # Define paths
    sample_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}_sorted.bam")
    spikein_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "spikein", f"{Configuration.file_to_process}_spikein_sorted.bam")
    stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")

    # Ensure stats directory exists
    os.makedirs(stats_dir, exist_ok=True)

    # Define output files
    flagstat_file = os.path.join(stats_dir, f"{Configuration.file_to_process}.flagstat.txt")
    samtools_stats_file = os.path.join(stats_dir, f"{Configuration.file_to_process}.samtools_stats.txt")
    idxstats_file = os.path.join(stats_dir, f"{Configuration.file_to_process}.idxstats.txt")

    # Check if BAM file exists
    if not os.path.exists(sample_bam):
        logging.error(f"BAM file {sample_bam} not found. Ensure alignment was completed before running stats.")
        raise Exception("BAM file not found")

    logging.info(f"Generating alignment stats for {Configuration.file_to_process}")

    # Run samtools flagstat
    with open(flagstat_file, "w") as outfile:
        subprocess.run(["samtools", "flagstat", sample_bam], stdout=outfile)
    logging.info(f"Samtools flagstat completed. Output saved to {flagstat_file}")

    # Run samtools stats
    with open(samtools_stats_file, "w") as outfile:
        subprocess.run(["samtools", "stats", sample_bam], stdout=outfile)
    logging.info(f"Samtools stats completed. Output saved to {samtools_stats_file}")

    # Run samtools idxstats
    with open(idxstats_file, "w") as outfile:
        subprocess.run(["samtools", "idxstats", sample_bam], stdout=outfile)
    logging.info(f"Samtools idxstats completed. Output saved to {idxstats_file}")

    logging.info("Finished generating alignment stats.")

    # If spike-in alignment exists, run samtools stats on spike-in BAM
    if Configuration.spike_in == True:
        spikein_stats_file = os.path.join(stats_dir, f"{Configuration.file_to_process}.spikein_stats.txt")
        with open(spikein_stats_file, "w") as outfile:
            subprocess.run(["samtools", "stats", spikein_bam], stdout=outfile)
        logging.info(f"Samtools stats for spike-in BAM saved to {spikein_stats_file}.")
    else:
        logging.info("[info] Skipping spike-in stats as spike-in alignment was not performed.")

