########################################
# script for merging fastq files, if needed



########################################

import os
import logging
import subprocess
from collections import defaultdict

def safe_symlink(src, dest):
    """Create a symlink, replacing existing links/files if necessary."""
    if os.path.islink(dest) or os.path.exists(dest):
        logging.warning(f"⚠️  Removing existing file/link at {dest}")
        os.remove(dest)
    logging.info(f"🔗 Creating symlink: {dest} → {src}")
    os.symlink(src, dest)

def merge_and_link_fastq_files(Configuration):
    """
    - Merges FASTQ files for samples with multiple runs/libraries.
    - Creates symbolic links in Processing_dir/sampleX/, whether merged or not.
    """

    # Define input directory for the current sample
    sample_raw_dir = os.path.join(Configuration.RAW_input_dir, Configuration.file_to_process)
    RAW_files = [f for f in os.listdir(sample_raw_dir) if f.endswith(".gz")]

    print("🔍 RAW Files Found:", RAW_files)

    # Group files by sample prefix
    sample_dict = defaultdict(lambda: {"R1": [], "R2": []})
    print("📂 Sample Dictionary Before Processing:", sample_dict)

    for file in RAW_files:
        sample_name = Configuration.file_to_process  # assuming one sample at a time
        if "_1.fq.gz" in file:
            sample_dict[sample_name]["R1"].append(file)
        elif "_2.fq.gz" in file:
            sample_dict[sample_name]["R2"].append(file)

    print("📂 Sample Dictionary After Processing:", sample_dict)

    # Define output directories
    sample_merged_dir = os.path.join(Configuration.Merged_dir, Configuration.file_to_process)
    sample_processing_dir = os.path.join(Configuration.Processing_dir, Configuration.file_to_process)

    os.makedirs(sample_merged_dir, exist_ok=True)
    os.makedirs(sample_processing_dir, exist_ok=True)

    for sample, reads in sample_dict.items():
        merged_R1 = os.path.join(sample_merged_dir, f"{sample}_R1.fastq.gz")
        merged_R2 = os.path.join(sample_merged_dir, f"{sample}_R2.fastq.gz")

        raw_R1_files = [os.path.join(sample_raw_dir, f) for f in reads["R1"]]
        raw_R2_files = [os.path.join(sample_raw_dir, f) for f in reads["R2"]]

        processing_R1 = os.path.join(sample_processing_dir, f"{sample}_R1.fastq.gz")
        processing_R2 = os.path.join(sample_processing_dir, f"{sample}_R2.fastq.gz")

        # --- Case 1: multiple R1/R2 files → merge them ---
        if len(reads["R1"]) > 1 or len(reads["R2"]) > 1:
            logging.info(f"📦 Multiple FASTQ files detected for {sample}. Merging...")

            with open(merged_R1, "wb") as outfile:
                subprocess.run(["cat"] + raw_R1_files, stdout=outfile, check=True)

            with open(merged_R2, "wb") as outfile:
                subprocess.run(["cat"] + raw_R2_files, stdout=outfile, check=True)

            safe_symlink(merged_R1, processing_R1)
            safe_symlink(merged_R2, processing_R2)

        # --- Case 2: single R1/R2 file → link directly ---
        else:
            logging.info(f"📎 Single FASTQ pair for {sample}. Linking directly.")
            safe_symlink(raw_R1_files[0], processing_R1)
            safe_symlink(raw_R2_files[0], processing_R2)

    logging.info("✅ Finished merging and linking FASTQ files.")

