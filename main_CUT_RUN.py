########################################
# main script that calls all functions

# naming system:
# folder name needs to be sampleName
# within the folder name needs to be id+lane_R1.fastq.gz and R2


# designed to work specifically for our library prepration methods
# only usable with paired end reads

########################################

from configuration import Config
import os
import glob
from random import random
from time import sleep
import argparse
import logging
from steps import merge, qc_raw_trimmed, trim_galore, alignment, filter_deduplication, bam_to_bedgraph, peak_calling, peak_calling_tester, peak_qc, FRiP


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Wrapper function for all steps of CUTnRUN analysis')

    parser.add_argument("-i",'--input', dest='infile', action='store', required=False,
                        help='input folder to force. Will overwrite all ouputs')
    parser.add_argument("-s",'--steps', dest='step', action='append', required=False,
                        help='chose steps instead of running everything')

    # parse arguments
    args = parser.parse_args()

    # set up configuration object for all steps. this sets up logging as well
    Configuration = Config()

    if args.infile == None:
        all_raws_present = [os.path.basename(x) for x in glob.glob(Configuration.RAW_input_dir)]

        all_processed = [os.path.basename(x) for x in glob.glob(Configuration.cleaned_alignments_dir)]
        # chose the first one of the ones that are still not processed and run 
        for i in all_raws_present:
            if i not in all_processed:
                os.makedirs(os.path.join(Configuration.cleaned_alignments_dir,i),exist_ok=True)
                Configuration.file_to_process = i
                break

        if Configuration.file_to_process == None:
            logging.error("There were no new files to process")
            raise Exception
        
    else:
        Configuration.file_to_process = args.infile
        os.makedirs(os.path.join(Configuration.cleaned_alignments_dir,Configuration.file_to_process),exist_ok=True)
    
    logging.info(f"This script will run the file : {Configuration.file_to_process}")

    if args.step == None:
        # call merging
        merge.merge_and_link_fastq_files(Configuration)

        # call qc
        qc_raw_trimmed.run_post_merge_qc(Configuration)
        
        # call trimming
        trim_galore.run_trim_galore(Configuration)
        
        # call post trim qc
        qc_raw_trimmed.run_post_trim_qc(Configuration)
        
        # call alignment
        alignment.run_bowtie2_alignment(Configuration)
        
        # call collect alignment stats
        alignment.generate_alignment_stats(Configuration)
        
        # filter aligned reads
        filter_deduplication.filter_reads(Configuration)
        filter_deduplication.spike_in_filter(Configuration)
        
        # mark duplicates (remove if IgG)
        filter_deduplication.mark_and_remove_duplicates(Configuration)
        
        # calculate normalisation
        filter_deduplication.calculate_normalisation(Configuration)
        
        # convert to bedgraph and bigwig
        bam_to_bedgraph.bam_to_bigwig(Configuration)

        # convert to bed
        bam_to_bedgraph.bam_to_bed(Configuration)
        
        # peak calling
        peak_calling.run_peak_calling(Configuration)
        
        # peak calling tester
        peak_calling_tester.run_sicer(Configuration)
        
        # peak qc
        peak_qc.FRiP(Configuration)
        
        

    else:
        if "merge" in args.step:
            merge.merge_and_link_fastq_files(Configuration)
        if "qc_raw" in args.step:
            qc_raw_trimmed.run_post_merge_qc(Configuration)
        if "trim_fastp" in args.step:
            trim_galore.run_fastp_trim(Configuration)
        if "qc_trim" in args.step:
            qc_raw_trimmed.run_post_trim_qc(Configuration)
        if "align" in args.step:
            alignment.run_bowtie2_alignment(Configuration)
        if "align_metrics" in args.step:
            alignment.generate_alignment_stats(Configuration)
        if "filter" in args.step:
            filter_deduplication.filter_reads(Configuration)
            filter_deduplication.spike_in_filter(Configuration)
        if "mark_dedup" in args.step:
            filter_deduplication.mark_and_remove_duplicates(Configuration)
        if "normalise" in args.step:
            filter_deduplication.calculate_normalisation(Configuration)
        if "convert" in args.step:
            bam_to_bedgraph.bam_to_bed(Configuration)
        if "call_peaks" in args.step:
            peak_calling.run_macs2_peak_calling(Configuration)
        if "test_SEACR" in args.step:
            peak_calling.run_peak_calling_v2(Configuration)
        if "test_peaks" in args.step:
            peak_calling_tester.run_sicer_peak_calling(Configuration)
        if "FRiP" in args.step:
            FRiP.FRiP(Configuration)
        if "peak_qc" in args.step:
            peak_qc.FRiP(Configuration)
