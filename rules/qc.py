# QUALITY CONTROL

# trim and screen ---------------------------------------------------------------------------------

rule fastp:
	input:
		r1 = "data/raw/{sample}_R1.fastq.gz",
		r2 = "data/raw/{sample}_R2.fastq.gz"
	output:
		r1 = "data/fastp/{sample}_R1.fastq.gz",
		r2 = "data/fastp/{sample}_R2.fastq.gz"
	conda:
		"../envs/fastp.yaml"
	log:
		"data/logs/{sample}.fastp.json"
	threads: 8
	shell:
		"fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
		"--detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

rule fastqc:
	input:
		"data/raw/{read}.fastq.gz"
	output:
		"data/fastqc/{read}_fastqc.html"
	conda:
		"../envs/fastqc.yaml"
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
		"../envs/fastq_screen.yaml"
	log:
		"data/logs/fastq_screen_{read}.txt"
	threads: 8
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen "
		"--conf {input.config} --force {input.fastq} > {log} 2>&1"

# library complexity ------------------------------------------------------------------------------

rule preseq:
	input:
		""
	output:
		"data/preseq/estimates_{sample}.txt"
	conda:
		"../envs/preseq.yaml"
	log:
		"data/logs/preseq_{sample}.log"
	shell:
		"preseq c_curve -P -o {output} {input} > {log} 2>&1" 

rule preseq_lcextrap:
	input:
		""
	output:
		"data/preseq/lcextrap_{sample}.txt"
	conda:
		"../envs/preseq.yaml"
	log:
		"data/logs/preseq_{sample}.log"
	shell:
		"preseq lc_extrap -P -e 1000000000 -o {output} {input} > {log} 2>&1"

# fragment length ---------------------------------------------------------------------------------

rule fraglength:
	input:
		"data/shift/{sample}.shifted.filtered.markd.sorted.bam"
	output:
		"data/stats/{sample}.fraglen.txt"
	conda:
		"../envs/bowtie2.yaml"
	shell:
		"samtools view {input} | cut -f9 | awk '$1 > 0 && $1 < 1000' | sort -S 4G | uniq -c | sort -b -k2,2n | awk '{{print $2, $1}}' > {output}"

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
			temp_df = pd.read_csv(i, sep = " ", index_col = 0, names = sample)
			dfs.append(temp_df)
		df = pd.concat(dfs, axis = 1)
		fraglen = df.plot()
		fraglen.update_layout( 
			title='Fragment Length Distribution', 
			xaxis_title='Fragment Length (bp)', 
			yaxis_title='Counts', 
			legend_title_text='Samples')
		fraglen.write_html(str(output))

# FRiP --------------------------------------------------------------------------------------------

# more like fraction of reads in consensus peaks

rule FRiP_count:
	input:
		consensus = "data/macs2/consensus_peaks.bed",
		sample = "data/shift/{sample}.shifted.filtered.markd.sorted.bam"
	output:
		"data/stats/{sample}.frip.txt"
	conda:
		"../envs/bedtools.yaml"
	shell:
		"""
		all_reads=$(samtools view -c {input.sample})
		rip=$(bedtools intersect -u -a {input.sample} -b {input.consensus} -ubam | samtools view -c)
		echo -e "{wildcards.sample}\n$all_reads\n$rip" > {output}
		"""

rule FRiP:
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

# multiqc -----------------------------------------------------------------------------------------

rule multiqc:
	input:
		directory("data/deseq2")
	output:
		"data/multiqc/multiqc_report.html"
	conda:
		"../envs/multiqc.yaml"
	shell:
		"multiqc -f data/ -o data/multiqc"