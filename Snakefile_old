#ATAC_seq2 pipeline rebooted
#After creating somehing
#Garth Kong and Rowan Callahan
#Maxson Lab 09/09/2020

#lets locate and find our configuration file
configfile: "config.yaml"
adapters_file = config["adapters_file"]
reference_align_genome = config["alignment_reference"] #The reference to align to
blacklist_file = config["blacklist_file"] if ( "blacklist_file" in config.keys() ) else "/home/groups/MaxsonLab/indices/GRch38/hg38.blacklist.v2.bed"
macs_genome_size = config["macs_genome_size"]
genome_name = config["genome_name"]
metadata_file = config["metadata_file"]

##now lets find the samples and create a list of them
SAMPLES, = glob_wildcards(raw_sample_folder + prefix + "{sample}" + readreverse_postfix)

import re
import pandas as pd

#we also want to make some files which are based on condition and not sample
#we now need to input our sample file to allow us to map condition to sample
config = pd.read_csv(metadata_file, delimiter='\t')
condition_to_samples_dictionary = {}
for sample in SAMPLES: 

    sample_id = sample 
    
    condition = config[config['SampleID'] == sample_id]['Condition'].values

    if str(*condition) in condition_to_samples_dictionary.keys():
        condition_to_samples_dictionary[str(*condition)].append(sample)
    else:
        condition_to_samples_dictionary[str(*condition)] = []
        condition_to_samples_dictionary[str(*condition)].append(sample)

SAMPLE_CONDITIONS = list(condition_to_samples_dictionary.keys()) 
print(SAMPLE_CONDITIONS)


rule all:
    input:
        expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_smooth.bw", merged_sample=MERGED_SAMPLES),
        expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_rough.bw", merged_sample=MERGED_SAMPLES),
        #This is the place holder, you cannot have multiple empty lines here otherwise you will run into errors
        expand(sample_work_path + "/fully_filtered/{merged_sample}_count_peaks_done.txt", merged_sample=MERGED_SAMPLES),
        expand(sample_work_path + "/bamfiles/reads_catalog_intervals_nodownsample.bed_{sample_condition}.bed", sample_condition=SAMPLE_CONDITIONS),
        expand(sample_work_path + "/fully_filtered/{sample_condition}_merged.bam", sample_condition = SAMPLE_CONDITIONS),
        sample_work_path + "/fully_filtered/all_read_catalog_counts.bed",
        sample_work_path + "/bamfiles/reads_catalog_intervals_nodownsample.bed",
        sample_work_path + "/bamfiles/intersected_bedfile_nodownsample.bed",
        sample_work_path + "/fully_filtered/intersected_bedfile_nodownsample_peakinfo.txt",
        sample_work_path + "/bamfiles/intersected_catalog.bed",
        sample_work_path + "/fully_filtered/all_read_catalog_nodownsample_counts.bed",
        sample_work_path + "/fully_filtered/correllation_plot.png",

        #sample_work_path + "/multiqc/done.txt",


include: "rules/filter_shift_downsample.smk"
include: "rules/peak_catalog_no_downsample.smk"
include: "rules/peak_catalog.smk"
#include: "rules/raw_bigwig.smk"
include: "rules/non_downsampled_bigwig.smk"
include: "rules/final_qc_metrics.smk"
include: "rules/merged_bigwig.smk"
include: "rules/quality_and_align.smk"
