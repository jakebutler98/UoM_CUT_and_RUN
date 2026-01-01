########################################
# utility that contains a class that gets the configuration



########################################

import logging
import datetime

class Config:
    """
    Class containing the parameters 
    """
    
    def __init__(self):

        self.RAW_input_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/raw/test_full_run"
        self.cleaned_alignments_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/clean/clean_alignments"
        self.Merged_dir= "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/merged"
        self.Processing_dir= "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/processing"
        self.RAW_QC_dir= "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/Raw_QC/"
        self.Trimmed_dir = "/mnt/iusers01/jw01/x25633jb/scratch/OA_CUTRUN_output/trimmed/"
        self.Trimmed_QC_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/trimmed_QC/"
        self.Alignment_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/Jake_Butler/OA_temp_storage/OA_CUTRUN_output/Alignments"
        self.Bowtie2_idx_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/common_files/data/external/Bowtie2_index_no_alt/GRCh38_no_alt"
        self.spikein_genome = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_Tag/common_files/data/external/spike_in/ecoli_K_12_MG1655_bowtie_index" ## note yet known - ASK ANA
        self.blacklist_bed = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/common_files/data/external/Blacklist/lists/hg38-blacklist.v2.bed"
        self.chrom_sizes = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/common_files/data/external/reference/Sequence/WholeGenomeFasta/chrom.sizes"
        self.MACS2_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/MACS2"
        self.SEACR_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/SEACR"
        self.peak_qc_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/peak_qc"
        self.SICER_dir = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/projects/CUT_RUN/analyses/CUT_RUN_processing/data/output/SICER"
        self.genome_chrom_sizes = "/mnt/jw01-aruk-home01/projects/oa_functional_genomics/common_files/data/external/hg38.chrom.sizes"
        
        self.spike_in = True        
        self.mapping_quality_threshold = "20"
        self.num_cores = "8"
        self.remove_mitochondrial = True
        

        self.bowtie2_index = "/mnt/jw01-aruk-home01/projects/shared_resources/sequencing/data/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
        self.genome_fasta = "/mnt/jw01-aruk-home01/projects/shared_resources/sequencing/data/Homo_sapiens/UCSC/hg38/Sequence/hg38/fasta/genome.fa"
        self.picard = "/mnt/jw01-aruk-home01/projects/functional_genomics/bin/picard/picard.jar"
        self.logs_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/ATAC_ChIP_pipeline/logs"

        self._init_logging()

        self.file_to_process = None
        self.analysis_type = None
        self.input_background = None
        
    def _init_logging(self):
        cur_date = datetime.datetime.now()
        
        logging.basicConfig(
            level=logging.INFO,
            format="%(levelname)s - %(message)s",
            handlers=[
                # logging.FileHandler("{0}/{1}.log".format(self.logs_dir, f"{cur_date.year}-{cur_date.month}-{cur_date.day}_{cur_date.hour}.{cur_date.minute}.{cur_date.second}"), mode="a"),
                logging.StreamHandler()]) 
