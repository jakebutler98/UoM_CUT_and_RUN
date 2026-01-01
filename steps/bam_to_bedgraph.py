# Converts BAM to BedGraph and BigWig using bamCoverage with spike-in normalization

import os
import logging
import subprocess

def bam_to_bigwig(Configuration):
    """
    Converts BAM to BedGraph and BigWig for both dupMarked and rmDup BAMs
    using deepTools' bamCoverage with spike-in normalization.
    Assumes BAM is sorted and indexed.
    """

    def run_bam_to_bigwig(file_type):
        input_bam = os.path.join(
            Configuration.Alignment_dir,
            Configuration.file_to_process,
            f"{Configuration.file_to_process}.{file_type}.bam"
        )
        bedgraph_file = os.path.join(
            Configuration.Alignment_dir,
            Configuration.file_to_process,
            f"{Configuration.file_to_process}.{file_type}.normalised.bedgraph"
        )
        bigwig_file = os.path.join(
            Configuration.Alignment_dir,
            Configuration.file_to_process,
            f"{Configuration.file_to_process}.{file_type}.normalised.bw"
        )

        spikein_stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")
        norm_file = os.path.join(spikein_stats_dir, f"{Configuration.file_to_process}.normalisation.txt")
        
        if not os.path.exists(norm_file):
            logging.error(f"Normalization factor file not found: {norm_file}")
            raise Exception("Normalization factor missing")

        if not os.path.exists(Configuration.chrom_sizes):
            logging.error(f"Chrom sizes file not found: {Configuration.chrom_sizes}")
            raise Exception("Missing chrom sizes file")

        # Read scale factor from last line of normalisation file
        with open(norm_file) as f:
            scale_factor = f.readlines()[-1].split(":")[1].strip()

        logging.info(f"Generating scaled BedGraph for {file_type}")

        # BedGraph
        bedgraph_cmd = [
            "bamCoverage",
            "--bam", input_bam,
            "--outFileName", bedgraph_file,
            "--outFileFormat", "bedgraph",
            "--normalizeUsing", "None",
            "--scaleFactor", scale_factor,
            "--numberOfProcessors", str(Configuration.num_cores)
        ]
        subprocess.run(bedgraph_cmd, check=True)

        logging.info(f"Scaled BedGraph created: {bedgraph_file}")

        # BigWig
        logging.info(f"Generating scaled BigWig for {file_type}")
        bigwig_cmd = [
            "bamCoverage",
            "--bam", input_bam,
            "--outFileName", bigwig_file,
            "--outFileFormat", "bigwig",
            "--normalizeUsing", "None",
            "--scaleFactor", scale_factor,
            "--numberOfProcessors", str(Configuration.num_cores)
        ]
        subprocess.run(bigwig_cmd, check=True)

        logging.info(f"Scaled BigWig created: {bigwig_file}")

        return bedgraph_file, bigwig_file

    # Process both rmDup and dupMarked BAMs
    bedgraph_rmdup, bigwig_rmdup = run_bam_to_bigwig("rmDup")
    bedgraph_dupmarked, bigwig_dupmarked = run_bam_to_bigwig("dupMarked")

    logging.info("BAM to BigWig conversion completed for both rmDup and dupMarked.")
    return {
        "rmDup": {"bedgraph": bedgraph_rmdup, "bigwig": bigwig_rmdup},
        "dupMarked": {"bedgraph": bedgraph_dupmarked, "bigwig": bigwig_dupmarked}
    }



def bam_to_bed(Configuration):
    """
    Converts BAM to BED using bedtools.
    Assumes BAM is sorted and indexed.
    """
    def run_command(command, log_message):
        subprocess.run(command, check=True)
        logging.info(log_message)

    def process_bam_to_bed(file_type):
        input_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.bam")
        temp_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.temp.bam")
        output_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.mapped.sorted.bam")
        clean_bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.mapped.sorted.fixmate.bam")
        bed_file = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.bed")
        clean_bed = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.bowtie2.bed")
        final_bed = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.bowtie2.fragments.bed")

        # Sort bam file (fixes pair-end read warning - heavily suggest this step so not to make 
        # log files MBs large)
        sort_cmd = [
            f"samtools view -u -f 0x2 -F 0x904 {input_bam} | samtools sort -n -@ 4 -m 1G -T {temp_bam} -o {output_bam} -"
        ]
        subprocess.run(sort_cmd, shell=True, check=True)
        # remove temp file
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
        logging.info(f"Sorted BAM created: {output_bam}")

        bam_to_bed_cmd = [
            f'bedtools bamtobed -bedpe -i {output_bam} | awk -v OFS="\\t" \'{{s=$2; e=$6; if(e<s){{t=s;s=e;e=t}} print $1,s,e}}\' > {clean_bed}'
        ]
        subprocess.run(bam_to_bed_cmd, shell=True, check=True)

        # sort by chr and pos
        bed_sort_cmd = [
            f"sort -k1,1 -k2,2n {clean_bed} > {final_bed}"
        ]
        subprocess.run(bed_sort_cmd, shell=True, check=True)
        logging.info(f"Final BED file created: {final_bed}")

        ### CONVERSION TO BEDGRAPH NORMALIZED TO SPIKE-IN USING NORMALISATION FACTOR ###
        norm_file = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats", f"{Configuration.file_to_process}.normalisation.txt")
        with open(norm_file) as f:
            scale_factor = f.readlines()[-1].split(":")[1].strip()
        
        logging.info(f"Generating normalised bedgraph (as input for SEACR) for {file_type}")
        bedgraph_file = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.{file_type}.bowtie2.fragments.normalised.bedgraph")
        bedtools_genomecov_cmd = [
            "bedtools", "genomecov",
            "-bg", "-scale", scale_factor,
            "-i", final_bed, "-g", Configuration.chrom_sizes,
        ]
        with open(bedgraph_file, "w") as out:
            subprocess.run(bedtools_genomecov_cmd, stdout=out, check=True)
        logging.info(f"Normalised bedgraph created: {bedgraph_file}")


    # Process both dedup and dupMarked files
    final_bed_dedup = process_bam_to_bed("rmDup")
    final_bed_dupmarked = process_bam_to_bed("dupMarked")

    logging.info("BAM to BED conversion completed.")



    return final_bed_dedup, final_bed_dupmarked

# test function for running bedtools genomecov
def test_bedtools_genomecov(Configuration):
    """
    Test function to run bedtools genomecov to create a bedgraph from a BAM file.
    """
    bam = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dupMarked.bam")
    chrom_sizes = Configuration.chrom_sizes
    output_bedgraph = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, f"{Configuration.file_to_process}.dupMarked.test.bedgraph")

    spikein_stats_dir = os.path.join(Configuration.Alignment_dir, Configuration.file_to_process, "stats")
    norm_file = os.path.join(spikein_stats_dir, f"{Configuration.file_to_process}.normalisation.txt")
        
    if not os.path.exists(norm_file):
        logging.error(f"Normalization factor file not found: {norm_file}")
        raise Exception("Normalization factor missing")

    if not os.path.exists(Configuration.chrom_sizes):
        logging.error(f"Chrom sizes file not found: {Configuration.chrom_sizes}")
        raise Exception("Missing chrom sizes file")

    # Read scale factor from last line of normalisation file
    with open(norm_file) as f:
        scale_factor = f.readlines()[-1].split(":")[1].strip()


    cmd = [
        "bedtools", "genomecov",
        "-ibam", bam,
        "-bg",
        "-g", chrom_sizes,
        "-scale", scale_factor
    ]
    with open(output_bedgraph, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    logging.info(f"Created bedgraph using bedtools genomecov: {output_bedgraph}")
    return output_bedgraph

