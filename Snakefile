""" Analyze your PE ATAC-Seq data """

import os
import sys
import glob
import pandas as pd
import plotly as plt
import plotly.graph_objects as go
from snakemake.utils import min_version
min_version("5.11")


SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")

def message(msg):
	sys.stderr.write("|--- " + msg + "\n")

for i in SAMPLES:
	message("Processing " + i)

configfile: "config.yaml"

if not os.path.isdir("data/stats"):
	os.mkdir("data/stats")

all_samples = glob.glob("data/raw/*.fastq.gz")
all_reads = [os.path.basename(i).split(".")[0] for i in all_samples]

localrules: fraglength_plot, FRiP, counts_table, multiqc

rule all:
	input:
		# quality control
		expand("data/fastp/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = all_reads),
		expand("data/fastq_screen/{sample}_{read}_screen.txt", sample = SAMPLES, read = ["R1", "R2"]),
		# expand("data/preseq/lcextrap_{sample}.txt", sample = SAMPLES),
		"data/fraglen.html",
		"data/frip.html",
		# read alignment
		expand("data/shift/{sample}.shifted.filtered.markd.sorted.bam", sample = SAMPLES),
		expand("data/bigwig/{sample}.bw", sample = SAMPLES),
		# peak calling
		expand("data/macs2/{sample}/{sample}_peaks.broadPeak", sample = SAMPLES),
		"data/macs2/consensus_peaks.bed",
		"data/counts/counts_table.txt",
		# differential
		directory("data/deseq2"),
		directory("data/diffbind"),
		"data/multiqc/multiqc_report.html",
		# HOMER
		directory("data/homer")

include: "rules/functions.py"
include: "rules/qc.py"
include: "rules/align.py"
include: "rules/peaks.py"
include: "rules/differential.py"
