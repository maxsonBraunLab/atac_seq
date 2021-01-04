#ATAC_seq2 pipeline rebooted

#Garth Kong and Rowan Callahan
#Maxson Lab 09/09/2020

import os
import pandas as pd
import plotly as plt
import plotly.graph_objects as go

configfile: "config.yaml"

SAMPLE, dir = glob_wildcards("samples/raw/{sample}_{dir,R1|R2}.fastq.gz")
SAMPLE = sorted(set(SAMPLE))
DIR = sorted(set(dir)) # R1, R2
# sample will glob cond+replicate information like 'MOLM24D_1' but not R1 or R2.

for i in SAMPLE:
	print("--- samples to process: {}".format(i))

essential_report = []
if config["gen_report"] == True: essential_report.append("data/essential_report.html"); print("Generating essential report")

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --profile ./slurm --cluster-config cluster.yaml

DB, CELL_LINE, FILE, = glob_wildcards("/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/{file}.bed.gz")

localrules: fragment_length_plot, frip_plot, counts_table, align_stats, essential_report

rule all:
	input:
		# trim, align reads, get fragment length distribution
		"samples/sample_stats/align_stats.txt",
		expand("samples/align/sorted/{sample}_sorted.bam", sample = SAMPLE),
		"data/fragment_length_dist.html",
		# fastqc and fastq_screen
		expand("samples/fastqc/{sample}_{dir}_paired_fastqc.{ext}", 
			sample = SAMPLE, dir = DIR, ext = ["html", "zip"]),
		expand("samples/fastq_screen/{sample}/{sample}_{dir}_paired_screen.{ext}", 
			sample = SAMPLE, dir = DIR, ext = ["png", "txt", "html"]),
		# filtered bamfiles, bigwigs, frip, counts table, consensus peaks
		"samples/macs/consensus_stats.txt",
		"data/frip.html",
		"data/counts_table.txt",
		expand("samples/bamfiles/filtered/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam", sample = SAMPLE),
		expand("data/bigwigs/{sample}_tracks.bw", sample = SAMPLE),
		# deseq2 for PCA and differential peak calc
		directory("data/de"),
		"data/de/norm_counts.txt",
		"data/de/log_norm_counts.txt",
		"data/de/sample_PCA.png",
		"data/de/de_stats.txt",
		essential_report,
		# chip screen - intersect consensus peaks with public chip data
		expand("data/chip_screen/{db}/{cell_line}/consensus_peaks_{file}.{ext}",
				zip, db = DB, cell_line = CELL_LINE, file = FILE, ext = ['bed.gz', 'out'] )


include: "rules/quality_and_align.smk"
include: "rules/filter_shift.smk"
include: "rules/peak_catalog_no_downsample.smk"
include: "rules/deseq2_report.smk"