import os
import logging
import subprocess

# build fastp function
def run_fastp_trim(Configuration):
    """
    Runs fastp on symlinked FASTQ files in the processing directory.
    Fastp will perform adapter trimming and quality filtering.
    """
    # Define input directory for the processing data (symlinked files)
    processing_sample_dir = os.path.join(Configuration.Processing_dir, Configuration.file_to_process)
    
    # Check if processing directory exists, if not, log and exit
    if not os.path.exists(processing_sample_dir):
        logging.error(f"Processing directory {processing_sample_dir} not found. Please run the merging, linking, and QC process first.")
        raise Exception("Processing directory not found")

    # Get the symlinked files in the processing directory (R1 and R2)
    processing_RAW_files = os.listdir(processing_sample_dir)
    processing_RAW_files = [f for f in processing_RAW_files if f.endswith(".gz")]
    
    # Ensure that the files are paired (i.e., R1 and R2 exist for each sample)
    if len(processing_RAW_files) % 2 != 0:
        logging.error(f"Unpaired files detected in {processing_sample_dir}. Ensure R1 and R2 files are present.")
        raise Exception("Unpaired files detected")

    # Create output directory for fastp results
    fastp_output_dir = os.path.join(Configuration.Trimmed_dir, Configuration.file_to_process)
    os.makedirs(fastp_output_dir, exist_ok=True)

    # Log the start of fastp
    logging.info(f"Running fastp on symlinked data for {Configuration.file_to_process}")

    # Run fastp on each pair of files
    for file in processing_RAW_files:
        if file.endswith("_R1.fastq.gz"):
            r1_file = os.path.join(processing_sample_dir, file)
            r2_file = os.path.join(processing_sample_dir, file.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            
            # Run fastp
            output_prefix = os.path.join(fastp_output_dir, os.path.basename(file).replace("_R1.fastq.gz", ""))
            subprocess.run([
                "fastp",
                "-i", r1_file,
                "-I", r2_file,
                "-o", output_prefix + "_R1.fastq.gz",
                "-O", output_prefix + "_R2.fastq.gz",
                "--length_required", "15",
                "--detect_adapter_for_pe",
                "--correction",
                "--trim_poly_g",
                "--overrepresentation_analysis",
                "--html", "--json", 
                "--thread", "4"
            ])
            
            # Logging the result
            logging.info(f"fastp completed for {file} and its pair. Output saved to {fastp_output_dir}")
            
    logging.info("Finished running fastp on all paired files.")
    