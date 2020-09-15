#ATAC_seq2 pipeline rebooted

#Garth Kong and Rowan Callahan
#Maxson Lab 09/09/2020

import os
import pandas as pd
import plotly as plt

configfile: "config.yaml"

SAMPLE, dir = glob_wildcards("samples/raw/{sample}_{dir,R1|R2}.fastq.gz")
SAMPLE = sorted(set(SAMPLE))
DIR = sorted(set(dir)) # R1, R2
# sample will glob cond+replicate information like 'MOLM24D_1' but not R1 or R2.

for i in SAMPLE:
	print("--- samples to process: {}".format(i))
# Comment out this if making dag or rulegraph.

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.N} -o {cluster.o} -e {cluster.e} -t {cluster.t} -J {cluster.J} -c {threads} --mem={cluster.mem}" -s Snakefile

localrules: fragment_length_plot

rule all:
	input:
		# trim, align reads, get fragment length distribution
		expand("samples/align/sorted/{sample}_sorted.bam", sample = SAMPLE),
		"samples/align/fragment_length/fragment_length_dist.html",
		# fastqc and fastq_screen
		expand("samples/fastqc/{sample}_{dir}_paired_fastqc.{ext}", 
			sample = SAMPLE, dir = DIR, ext = ["html", "zip"]),
		expand("samples/fastq_screen/{sample}/{sample}_{dir}_paired_screen.{ext}", 
			sample = SAMPLE, dir = DIR, ext = ["png", "txt", "html"]),
        #filtered bamfiles 
        expand("samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam", sample = SAMPLE),


include: "rules/quality_and_align.smk"
include: "rules/filter_shift_downsample.smk"
