#ATAC_seq pipeline
#Rowan Callahan
#Maxson Lab 02/14/2020
#
#
#
#
#


configfile: "config.yaml"

# add prefix and then read 1 postfix read 2 postfix and then can have those in the config file
#prefix = config["prefix"] you can enter this in the config file or just in the snakemake file if you want
#read1_postfix = config["r1_postfix"]
#read2_postfix = config["r2_postfix"]
#sample_path = config["sample_path"]

prefix = config["prefix"] #"LIB200212TB_"
readforward_postfix = config["readforward_postfix"] #"R1_001.fastq.gz"
readreverse_postfix = config["readreverse_postfix"] #"R2_001.fastq.gz"
sample_work_path = config["sample_work_path"] #"../atac_run_20200220" #"../samples"
raw_sample_folder = config["raw_sample_folder"] #"/home/groups/MaxsonLab/input-data/CUT_RUNTAG/2020/LIB200212TB/200218_A01058_0021_AHKCHWDRXX/LIB200212TB/" #"../samples/raw/" #where the raw samples are kept
adapters_file = config["adapters_file"]
reference_align_genome = config["alignment_reference"] #"/home/groups/MaxsonLab/indices/GRch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
LANES = config["LANES"] #['_L001_','_L002_']

SAMPLES, = glob_wildcards(raw_sample_folder + prefix + "{sample}" + readreverse_postfix)
#MERGED_SAMPLES, = glob_wildcards(raw_sample_folder + prefix + "{merged_sample}"+ LANES[0] + readreverse_postfix)
MERGED_SAMPLES = SAMPLES

print(SAMPLES)
print(MERGED_SAMPLES)

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_smooth.bw", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_rough.bw", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",merged_sample=MERGED_SAMPLES),
        expand(sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv", merged_sample=MERGED_SAMPLES),

include: "rules/begin_qc.smk"
include: "rules/postalign_qc.smk"


