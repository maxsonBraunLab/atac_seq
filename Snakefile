#ATAC_seq pipeline
#Rowan Callahan
#Maxson Lab 02/14/2020

#lets locate and find our configuration file
configfile: "config.yaml"

#now lets load all of our configurations from the config file into variables in our script
prefix = config["prefix"] 
readforward_postfix = config["readforward_postfix"] 
readreverse_postfix = config["readreverse_postfix"] 

sample_work_path = config["sample_work_path"] 
raw_sample_folder = config["raw_sample_folder"] #where the raw samples are kept

adapters_file = config["adapters_file"]
reference_align_genome = config["alignment_reference"] #The reference to align to

blacklist_file = config["blacklist_file"] if ( "blacklist_file" in config.keys() ) else "/home/groups/MaxsonLab/indices/GRch38/hg38.blacklist.v2.bed"
macs_genome_size = config["macs_genome_size"] if ( "macs_genome_size" in config.keys() ) else "2.7e9"
genome_name = config["genome_name"] if ( "genome_name" in config.keys() ) else "hg38"
fastq_screen_conf = "/home/groups/MaxsonLab/callahro/atac_seq/fastqscreen.conf"

metadata_file = config["metadata_file"] if ("metadata_file" in config.keys() ) else "/home/groups/MaxsonLab/callahro/atac_seq/MOLM_ATAC_metadata.txt"
remove_filepath_regex = config["remove_filepath_regex"] if ("remove_filepath_regex" in config.keys() ) else "-ATAC/.*"

LANES = config["LANES"] 

##now lets find the samples and create a list of them
SAMPLES, = glob_wildcards(raw_sample_folder + prefix + "{sample}" + readreverse_postfix)
MERGED_SAMPLES, = glob_wildcards(raw_sample_folder + prefix + "{merged_sample}"+ LANES[0] + readreverse_postfix)
print(raw_sample_folder + prefix + LANES[0] + readreverse_postfix)
#MERGED_SAMPLES = SAMPLES

print(SAMPLES)
print(MERGED_SAMPLES)

#now lets filter out the samples that we don't want and that we want to include from our sample list
include_samples_regex = config["include_samples_regex"] if ("include_samples_regex" in config.keys() ) else r'.*' 
exclude_samples_regex = config["exclude_samples_regex"] if ("exclude_samples_regex" in config.keys() ) else r' ' 

import re
import pandas as pd

#include the samples that you want to include and exclude the samples
print("----now filtering out the samples that we care about")
SAMPLES = [ sample for sample in SAMPLES if re.match(include_samples_regex, sample) and (not re.match(exclude_samples_regex, sample) ) ] 
print(SAMPLES)
MERGED_SAMPLES = [ sample for sample in MERGED_SAMPLES if re.match(include_samples_regex, sample) and (not re.match(exclude_samples_regex, sample) ) ] 
print(MERGED_SAMPLES)

#we also want to make some files which are based on condition and not sample
#we now need to input our sample file to allow us to map condition to sample
config = pd.read_csv(metadata_file, delimiter='\t')
condition_to_samples_dictionary = {}
for sample in MERGED_SAMPLES:

    sample_id = re.sub(remove_filepath_regex, "", str(sample))
    
    condition = config[config['SampleID'] == sample_id]['Condition'].values

    if str(*condition) in condition_to_samples_dictionary.keys():
        condition_to_samples_dictionary[str(*condition)].append(sample)
    else:
        condition_to_samples_dictionary[str(*condition)] = []
        condition_to_samples_dictionary[str(*condition)].append(sample)

SAMPLE_CONDITIONS = list(condition_to_samples_dictionary.keys()) 
print(SAMPLE_CONDITIONS)

#def samples_for_merged_bigwig(wildcards):
#    return([ str("/list/to/file/" + name + "/rest/file.txt") for name in condition_to_samples_dictionary[wildcards.merged_sample] ])


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        #expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/bamfiles/{merged_sample}_homerout/homerMotifs.all.motifs", merged_sample=MERGED_SAMPLES)
        #expand(sample_work_path + "/raw_bigwigs/{merged_sample}_raw.bw", merged_sample = MERGED_SAMPLES),
        #expand(sample_work_path + "/non_downsampled_bigwigs/{merged_sample}_non_downsampled.bw", merged_sample = MERGED_SAMPLES),
        #expand(sample_work_path + "/fastqc/{sample}_fastqc_done.txt", sample = MERGED_SAMPLES), 
        #expand(sample_work_path + "/fastqscreen/{sample}_fastqscreen_done.txt", sample = MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_information_done.txt", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_flagstat.txt", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_smooth.bw", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_rough.bw", merged_sample=MERGED_SAMPLES),
        #This is the place holder, you cannot have multiple empty lines here otherwise you will run into errors
        #expand(sample_work_path + "/fully_filtered/{merged_sample}_count_peaks_done.txt", merged_sample=MERGED_SAMPLES),
        #expand(sample_work_path + "/bamfiles/reads_catalog_intervals_nodownsample.bed_{sample_condition}.bed", sample_condition=SAMPLE_CONDITIONS),
        #expand(sample_work_path + "/fully_filtered/{sample_condition}_merged.bam", sample_condition = SAMPLE_CONDITIONS),
        #sample_work_path + "/fully_filtered/all_read_catalog_counts.bed",
        #sample_work_path + "/bamfiles/reads_catalog_intervals_nodownsample.bed",
        #sample_work_path + "/bamfiles/intersected_bedfile_nodownsample.bed",
        #sample_work_path + "/fully_filtered/all_read_catalog_nodownsample_counts.bed",
        #sample_work_path + "/fully_filtered/intersected_bedfile_nodownsample_peakinfo.txt",
        sample_work_path + "/bamfiles/intersected_catalog.bed",
        sample_work_path + "/fully_filtered/all_read_catalog_nodownsample_counts_ondownsamplepeaks.bed",

        #sample_work_path + "/multiqc/done.txt",


include: "rules/filter_shift_downsample.smk"
include: "rules/peak_catalog_no_downsample.smk"
include: "rules/peak_catalog.smk"
#include: "rules/raw_bigwig.smk"
include: "rules/non_downsampled_bigwig.smk"
include: "rules/final_qc_metrics.smk"
include: "rules/merged_bigwig.smk"
