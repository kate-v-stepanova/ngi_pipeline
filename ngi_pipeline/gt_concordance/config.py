import os
# ------------------------------------------------------------------------
# IMPORTANT !!!!
# the pathes must be configured exactly in the same format:
# when there is a var {sample} or {project} in path - it must remain there
# python script will use this variables and insert the values
# eg. ANALYSIS_PATH.format(project='P1234')
# if there is no {project} inside of the string, it will cause an Exception
# -------------------------------------------------------------------------


# Set path to where to find analysis for p-id
ANALYSIS_PATH = "/apus/v1/a2014205/nobackup/NGI/analysis_ready/ANALYSIS/{project}/piper_ngi"
BAM_FILE = os.path.join(ANALYSIS_PATH, "05_processed_alignments/{sample}.clean.dedup.recal.bam")

# Set path to where you find MAF_GT files
# WORK_DIR = "/apus/h1/funk_002/bin/GT_concordance" # file_path
WORK_DIR = "/proj/a2015179/nobackup/perlu/GT_concordance"
# where results will be written -> .gt file, .conc files, .vcf file
OUTPUT_DIR = "/apus/v1/a2014205/nobackup/NGI/analysis_ready/ANALYSIS/{project}/piper_ngi/08_misc/{sample}" # todo: must be PROJECT_DIR,
                                                                                                            # todo: SAMPLE_DIR
                                                                                                            # todo: output_dir = PROJECT_DIR + <project> + SAMPLE_DIR + <sample>
# OUTPUT_DIR = WORK_DIR
# reference file -> for java command
GATK_REF_FILE = "/proj/a2014205/piper_references/gatk_bundle/2.8/b37/human_g1k_v37.fasta"
# variants file -> for java command
GATK_VAR_FILE = "/proj/a2014205/piper_references/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf"
# in case if module load was not executed
GATK_HOME = "/sw/apps/bioinfo/GATK/3.4.0" # todo: setup virtual environment

# input files:
INTERVAL_FILE = os.path.join(WORK_DIR, "snps.interval_list")
SNPS_FILE = os.path.join(WORK_DIR, 'maf_snps.txt')


# concordance output file -> the result of check_concordance will be writen there
CONC_FILE = os.path.join(OUTPUT_DIR, '{sample}.conc')

# main concordance file -> aggregation of conc_files
MAIN_CONC_FILE = os.path.join(OUTPUT_DIR, '{project}.concordance')
