#Snakemake has a special function unpack that allows you to list functions as input rules
#These functoins can take in the wildcards object and allow you to create a list of files from a given wildcard input
#IN this case snakemake takes the output files that will be generated, one for each sample condition
#Then snakemake feeds that output {sample_condition} to the input of the sample_for_merged_bigwig function
#This function will generate a list of all of the files that our config file tells us is are member of this condition
def samples_for_merged_bigwig(wildcards):
    return([ str(sample_work_path + "/bamfiles/" + name + "_rmChrM_dedup_quality_shiftedReads_downSample.bam") for name in condition_to_samples_dictionary[wildcards.sample_condition] ]) 

def create_bigwig_sample_list(wildcards):
    return(" ".join([ str(sample_work_path + "/bamfiles/" + name + "_rmChrM_dedup_quality_shiftedReads_downSample.bam") for name in condition_to_samples_dictionary[wildcards.sample_condition] ]))

#This rule creates a merged bigwig file where the samples are merged by their condition
rule makeMergedBigwig:
    input:
        unpack(samples_for_merged_bigwig) 
    output:
        bamfile =  sample_work_path + "/fully_filtered/{sample_condition}_merged.bam",
        bigwig =  sample_work_path + "/fully_filtered/{sample_condition}_merged.bw",
    params:
        merge_input = lambda wildcards: create_bigwig_sample_list(wildcards) 
    conda:
        "../envs/deeptools.yaml"
    threads: 8 
    shell:
        """
        export PATH=$PATH:/opt/installed/samtools-1.6/bin/
        samtools merge -@ {threads} - {params.merge_input} > {output.bamfile}
        samtools index {output.bamfile}
        bamCoverage --bam {output.bamfile} -o {output.bigwig} --numberOfProcessors {threads} --skipNAs --binSize 5 --smoothLength 15  
        """

