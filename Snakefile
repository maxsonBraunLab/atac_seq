#ATAC_seq pipeline
#Rowan Callahan
#Maxson Lab 02/14/2020
#
#
#
#
#


configfile: "config.yaml"

SAMPLES, = glob_wildcards("../samples/raw/{sample}_1.fastq")

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        expand("../samples/samfiles/{sample}.sam", sample = SAMPLES),

include: "rules/begin_qc.smk"
##include: "rules/deseq.smk"



