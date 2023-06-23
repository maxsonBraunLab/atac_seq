""" Analyze your PE ATAC-Seq data """

import os
import sys
import glob
import pandas as pd
import plotly as plt
import plotly.graph_objects as go
from snakemake.utils import min_version
min_version("5.11")
if sys.version_info < (3, 6):
	sys.exit("Python version is less than 3.6. Your python version:", sys.version_info)

SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")

def message(msg):
	sys.stderr.write("|--- " + msg + "\n")

for i in SAMPLES:
	message("Processing " + i)

def defect_mode(wildcards, attempt):
	if attempt == 1:
		return ""
	elif attempt > 1:
		return "-D"

configfile: "config/config.yaml"

all_samples = glob.glob("data/raw/*.fastq.gz")
all_reads = [os.path.basename(i).split(".")[0] for i in all_samples]
all_conditions = set([ os.path.basename(i).split("_")[0] for i in all_samples ])

def get_tracks_by_condition(wildcards):
	samples_by_condition = [ i for i in SAMPLES if wildcards.condition in i ]
	mergebw_input = [ "data/bigwig/{}.bw".format(i) for i in samples_by_condition ]
	return(mergebw_input)

if config["ASSEMBLY"] == "hg38":
	GSIZE = 'hs'
elif config["ASSEMBLY"] == "mm10":
	GSIZE = 'mm'
else:
	sys.exit("ERROR: Only hg38 and mm10 are supported. Your assembly: " + config["ASSEMBLY"])

localrules: fraglength_plot, FRiP_plot, counts_table, multiqc

rule all:
	input:
		# quality control -------------------------------------------------------------------------
		expand("data/fastp/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = all_reads),
		expand("data/fastq_screen/{sample}_{read}_screen.txt", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/preseq/estimates_{sample}.txt", sample = SAMPLES),
		expand("data/preseq/lcextrap_{sample}", sample = SAMPLES),
		"data/multiqc/multiqc_report.html",
		"data/fraglen.html",
		"data/frip.html",
		# read alignment --------------------------------------------------------------------------
		expand("data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam", sample = SAMPLES),
		expand("data/bigwig/{sample}.bw", sample = SAMPLES),
		expand("data/mergebw/{conditions}.bw", conditions = all_conditions),
		# peak calling ----------------------------------------------------------------------------
		expand("data/macs2/{sample}_peaks.broadPeak", sample = SAMPLES),
		"data/counts/consensus_peaks.bed",
		"data/counts/counts_table.txt",
		# differential ----------------------------------------------------------------------------
		"data/deseq2",
		"data/diffbind",
		"data/homer"

# pre-processing ----------------------------------------------------------------------------------

rule fastp:
	input:
		r1 = "data/raw/{sample}_R1.fastq.gz",
		r2 = "data/raw/{sample}_R2.fastq.gz"
	output:
		r1 = "data/fastp/{sample}_R1.fastq.gz",
		r2 = "data/fastp/{sample}_R2.fastq.gz"
	conda:
		"envs/fastp.yaml"
	log:
		"data/logs/{sample}.fastp.json"
	threads: 8
	shell:
		"fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
		"--detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

rule fastqc:
	input:
		"data/fastp/{read}.fastq.gz"
	output:
		"data/fastqc/{read}_fastqc.html"
	conda:
		"envs/fastqc.yaml"
	log:
		"data/logs/fastqc_{read}.log"
	threads: 4
	shell:
		"fastqc -t {threads} --outdir data/fastqc {input} > {log} 2>&1"

rule fastq_screen:
	input:
		fastq = "data/fastp/{read}.fastq.gz",
		config = config["FASTQ_SCREEN_CONFIG"]
	output:
		"data/fastq_screen/{read}_screen.txt"
	conda:
		"envs/fastq_screen.yaml"
	log:
		"data/logs/fastq_screen_{read}.txt"
	threads: 8
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen "
		"--conf {input.config} --force {input.fastq} > {log} 2>&1"

rule bwa:
	input:
		r1 = rules.fastp.output.r1,
		r2 = rules.fastp.output.r2
	output:
		temp("data/bwa/{sample}.sorted.bam")
	conda:
		"envs/bwa.yaml"
	log:
		"data/logs/bwa_{sample}.log"
	threads: 8
	shell:
		"bwa mem -t {threads} {config[GENOME]} {input.r1} {input.r2} 2>{log} | samtools sort -@ {threads} > {output}"

rule filter:
	input:
		rules.bwa.output
	output:
		temp("data/filter/{sample}.filtered.sorted.bam")
	conda:
		"envs/bwa.yaml"
	threads: 4
	shell:
		"samtools view -@ {threads} -h -F 1804 -f 2 {input[0]} | "
		"grep -v chrM | samtools sort -@ {threads} > {output}"

# -F 1804 = filter away away 1804
# 	read paired (1), read unmapped (4),
# 	mate unmapped (8), not in primary alignment (100),
# 	read fails platform (200)
# 	no duplicate reads
# -f 2 = filter for paired-end reads

rule rmdup:
	input:
		rules.filter.output
	output:
		temp("data/rmdup/{sample}.rmdup.sorted.bam"),
		temp("data/rmdup/{sample}.rmdup.sorted.bam.bai")
	params:
		dup = "data/stats/{sample}.dup.txt"
	conda:
		"envs/sambamba.yaml"
	log:
		"data/logs/rmdup_{sample}.log"
	threads: 4 
	shell:
		"sambamba markdup -r -t {threads} --tmpdir=data/rmdup --io-buffer-size=512 {input} {output[0]} > {log} 2>&1"

rule banlist:
	input:
		bam = rules.rmdup.output[0],
		banlist = config["BANLIST"]
	output:
		"data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam",
		"data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam.bai"
	params:
		banlist = "data/stats/{sample}.banlist.txt"
	conda:
		"envs/bedtools.yaml"
	threads: 4
	shell:
		"bedtools intersect -v -ubam -abam {input.bam} -b {input.banlist} | samtools sort -@ {threads} > {output[0]}; samtools index {output[0]}"

rule bigwig:
	input:
		bam = rules.banlist.output[0],
		bai = rules.banlist.output[1],
	output:
		"data/bigwig/{sample}.bw"
	conda:
		"envs/deeptools.yaml"
	threads: 8
	shell:
		"bamCoverage -b {input[0]} -o {output} --binSize 10 --smoothLength 50 --normalizeUsing CPM -p {threads} "

rule mergebw:
	input:
		get_tracks_by_condition
	output:
		"data/mergebw/{condition}.bw"
	conda:
		"envs/mergebw.yaml"
	threads: 8
	shell:
		"bash scripts/mergebw.sh -c {config[CHROM_SIZES]} -o {output} {input}"

rule fraglength:
	input:
		"data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam"
	output:
		"data/stats/{sample}.fraglen.txt"
	conda:
		"envs/bowtie2.yaml"
	shell:
		"samtools view {input} | awk '$9>0 && $9 < 1000' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | awk -v OFS='\t' '{{print $2,$1}}' > {output}"

rule fraglength_plot:
	input:
		expand("data/stats/{sample}.fraglen.txt", sample = SAMPLES)
	output:
		"data/fraglen.html"
	run:
		pd.options.plotting.backend = "plotly"
		dfs = []
		for i in input:
			sample = [os.path.basename(i).split(".")[0]]
			temp_df = pd.read_csv(i, sep = "\t", index_col = 0, names = sample)
			dfs.append(temp_df)
		df = pd.concat(dfs, axis = 1)
		fraglen = df.plot()
		fraglen.update_layout( 
			title='Fragment Length Distribution', 
			xaxis_title='Fragment Length (bp)', 
			yaxis_title='Counts', 
			legend_title_text='Samples')
		fraglen.write_html(str(output))

# more like fraction of reads in consensus peaks
rule FRiP:
	input:
		consensus = "data/counts/consensus_peaks.bed",
		sample = "data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam"
	output:
		"data/stats/{sample}.frip.txt"
	conda:
		"envs/bedtools.yaml"
	shell:
		"""
		all_reads=$(samtools view -c {input.sample})
		rip=$(bedtools intersect -u -a {input.sample} -b {input.consensus} -ubam | wc -l)
		echo -e "{wildcards.sample}\n$all_reads\n$rip" > {output}
		"""

rule FRiP_plot:
	input:
		expand("data/stats/{sample}.frip.txt", sample = SAMPLES)
	output:
		"data/frip.html"
	run:
		pd.options.plotting.backend = "plotly"
		dfs = []
		for i in input:
			# sample = [os.path.basename(i).split(".")[0]]
			temp_df = pd.read_csv(i, sep = " ")
			dfs.append(temp_df)
		df = pd.concat(dfs, axis = 1)
		df = df.rename(index={0: 'total_reads', 1: 'reads_in_peaks'})
		df.loc['ratio'] = df.loc['reads_in_peaks'] / df.loc['total_reads']

		# plot graph. plot ratio as bottom as percent, and plot to max value of 1.
		fig = go.Figure(data=[
			go.Bar(name='inside_peaks', x=df.columns, y=df.loc['ratio'], marker_color='rgb(255,201,57)'),
			go.Bar(name='outside_peaks', x=df.columns, y= ([1] * df.shape[1]) - df.loc['ratio'], marker_color='rgb(0,39,118)')])
		# Change the bar mode
		fig.update_layout(barmode='stack', title='Fraction of Reads in Peaks by Sample', xaxis_tickfont_size=14,
			yaxis=dict(title='Fraction of reads in peaks', titlefont_size=16, tickfont_size=14),
			xaxis=dict(title='Samples'))
		fig.write_html(str(output))

rule preseq:
	input:
		rules.banlist.output[0]
	output:
		"data/preseq/estimates_{sample}.txt"
	conda:
		"envs/preseq.yaml"
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_{sample}.log"
	shell:
		"preseq c_curve -B {resources.defect_mode} -l 1000000000 -P -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
	input:
		rules.banlist.output[0]
	output:
		"data/preseq/lcextrap_{sample}"
	conda:
		"envs/preseq.yaml"
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_lcextrap_{sample}.log"
	shell:
		"preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -P -e 1000000000 -o {output} {input} > {log} 2>&1"

# peak calling ------------------------------------------------------------------------------------

rule macs2:
	input:
		bam = "data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam",
		bai = "data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam.bai"
	output:
		"data/macs2/{sample}_peaks.broadPeak"
	conda:
		"envs/macs2.yaml"
	params:
		GS = GSIZE
	log:
		"data/logs/macs2_{sample}.log"
	shell:
		"macs2 callpeak -t {input.bam} -n {wildcards.sample} "
		"--format BAMPE --gsize {params.GS} --tempdir data/macs2 "
		"--outdir data/macs2 --broad > {log} 2>&1"

rule consensus:
	input:
		expand("data/macs2/{sample}_peaks.broadPeak", sample = SAMPLES)
	output:
		"data/counts/consensus_peaks.bed"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bash scripts/consensus_peaks.sh {config[N_INTERSECTS]} {input}"

rule counts:
	input:
		sample = "data/banlist/{sample}.banlist.filtered.rmdup.sorted.bam",
		consensus = "data/counts/consensus_peaks.bed"
	output:
		"data/multicov/{sample}.txt"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools multicov -bams {input.sample} -bed {input.consensus} > {output}"

rule counts_table:
	input:
		expand("data/multicov/{sample}.txt", sample = sorted(SAMPLES))
	output:
		"data/counts/counts_table.txt"
	run:
		dfs = []
		for file in list(input):
			sample_name = os.path.basename(file).split('.')[0]
			dfs.append( pd.read_csv(file, sep = "\t", names = ["chr", "start", "end", sample_name]) )
		df = pd.concat(dfs, axis = 1)
		df = df.loc[:,~df.columns.duplicated()]
		df.insert( 3, "name", ["peak" + str(i) for i in df.index] )
		df.to_csv(str(output), header = True, index = False, sep = "\t")

# differential testing ----------------------------------------------------------------------------

rule deseq2:
	input:
		counts = "data/counts/counts_table.txt",
		metadata = config["DESEQ2_CONFIG"]
	params:
		padj_cutoff = config["padj_cutoff"]
	output:
		directory("data/deseq2"),
		norm_counts = "data/deseq2/norm_counts.txt",
		log_norm_counts = "data/deseq2/log_norm_counts.txt",
		pca = "data/deseq2/sample_PCA.png",
		stats = "data/deseq2/de_stats.txt",
		all_sig_intervals = "data/deseq2/all_sig_intervals.bed",
		contrast_combinations = "data/deseq2/contrast_combinations.txt"
	conda:
		"envs/deseq2.yaml"
	threads: 8
	log:
		out = "data/logs/deseq2.log"
	script:
		"scripts/deseq2.R"
# normalize by reads in peaks

rule diffbind:
	input:
		consensus_peaks = "data/counts/consensus_peaks.bed",
		metadata = config["DIFFBIND_CONFIG"]
	params:
		padj_cutoff = config["padj_cutoff"]
	output:
		directory("data/diffbind")
	conda:
		"envs/diffbind.yaml"
	threads: 8
	log:
		"data/logs/diffbind.log"
	script:
		"scripts/diffbind.R"
# normalize by entire sequencing depth

# homer motif analysis ----------------------------------------------------------------------------

rule homer:
	input:
		rules.deseq2.output.contrast_combinations
	output:
		directory("data/homer")
	params:
		genome = config["GENOME"]
	log:
		"data/logs/homer.log"
	conda:
		"envs/homer.yaml"
	shell:
		"bash scripts/homer.sh -i {input} -g {params.genome} -s 0 -c 8"
# this rule submits HOMER runs to SLURM if -s = 1. A run is each unique contrast
# combinations split by up and down peaks if DE peaks >= 10.

rule multiqc:
	input:
		expand("data/fastp/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = all_reads),
		expand("data/fastq_screen/{sample}_{read}_screen.txt", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/macs2/{sample}_peaks.broadPeak", sample = SAMPLES),
		expand("data/preseq/lcextrap_{sample}", sample = SAMPLES)
	output:
		"data/multiqc/multiqc_report.html"
	conda:
		"envs/multiqc.yaml"
	shell:
		"multiqc -f data/ -o data/multiqc --ignore data/homer"
