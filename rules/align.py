rule bwa:
	input:
		r1 = rules.fastp.output.r1,
		r2 = rules.fastp.output.r2
	output:
		temp("data/bwa/{sample}.sorted.bam")
	conda:
		"../envs/bwa.yaml"
	log:
		"data/logs/bwa_{sample}.log"
	threads: 8
	shell:
		"bwa mem -t {threads} {config[GENOME]} {input.r1} {input.r2} 2>{log} | samtools sort -@ {threads} -o {output} -"

rule mito:
	input:
		rules.bwa.output
	output:
		temp("data/mito/{sample}.mito.sorted.bam")
	params:
		mito = "data/stats/{sample}.mito.txt"
	conda:
		"../envs/bwa.yaml"
	threads: 4
	shell:
		"""
		with_mito=$(samtools view -@ {threads} -c {input})
		samtools view -h {input} | grep -v chrM | samtools sort -@ {threads} > {output}
		without_mito=$(samtools view -@ {threads} -c {output})
		echo -e "{wildcards.sample}\n$with_mito\n$without_mito" > {params.mito}
		"""

rule markdup:
	input:
		rules.mito.output
	output:
		temp("data/markd/{sample}.markd.sorted.bam"),
		temp("data/markd/{sample}.markd.sorted.bam.bai")
	params:
		dup = "data/stats/{sample}.dup.txt"
	conda:
		"../envs/sambamba.yaml"
	log:
		"data/logs/markdup_{sample}.log"
	threads: 4 
	shell:
		"""
		with_dup=$(samtools view -@ {threads} -c {input})
		sambamba markdup -r -t {threads} --tmpdir=data/markd --io-buffer-size=512 {input} {output[0]} > {log} 2>&1
		without_dup=$(samtools view -@ {threads} -c {output[0]})
		echo -e "{wildcards.sample}\n$with_dup\n$without_dup" > {params.dup}
		"""

rule quality:
	input:
		rules.markdup.output
	output:
		temp("data/quality/{sample}.filtered.markd.sorted.bam")
	params:
		quality = "data/stats/{sample}.quality.txt"
	conda:
		"../envs/bwa.yaml"
	threads: 4
	shell:
		"""
		all_reads=$(samtools view -@ {threads} -c {input[0]})
		samtools view -@ {threads} -h -b -F 1804 -f 2 {input[0]} -o {output}
		quality_reads=$(samtools view -@ {threads} -c {output})
		echo -e "{wildcards.sample}\n$all_reads\n$quality_reads" > {params.quality}
		"""
# filter away (1804) =
# 	read paired (1), read unmapped (4),
# 	mate unmapped (8), not in primary alignment (100),
# 	read fails platform (200)
# 	no duplicate reads

rule banlist:
	input:
		bam = rules.quality.output,
		banlist = config["BANLIST"]
	output:
		"data/banlist/{sample}.banlist.filtered.markd.sorted.bam",
		"data/banlist/{sample}.banlist.filtered.markd.sorted.bam.bai"
	params:
		banlist = "data/stats/{sample}.banlist.txt"
	conda:
		"../envs/bedtools.yaml"
	threads: 4
	shell:
		"""
		all_reads=$(samtools view -@ {threads} -c {input.bam})
		bedtools intersect -v -a {input.bam} -b {input.banlist} > {output[0]}
		samtools index {output[0]}
		good_reads=$(samtools view -@ {threads} -c {output[0]})
		echo -e "{wildcards.sample}\n$all_reads\n$good_reads" > {params.banlist}
		"""

rule shift:
	input:
		rules.banlist.output[0]
	output:
		temp("data/shift/{sample}.shifted.filtered.markd.bam")
	conda:
		"../envs/deeptools.yaml"
	threads: 4
	shell:
		"alignmentSieve -b {input} -o {output} -p {threads} --ATACshift"

rule sort:
	input:
		rules.shift.output
	output:
		"data/shift/{sample}.shifted.filtered.markd.sorted.bam"
	conda:
		"../envs/bwa.yaml"
	threads: 4
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

rule index:
	input:
		rules.sort.output
	output:
		"data/shift/{sample}.shifted.filtered.markd.sorted.bam.bai"
	conda:
		"../envs/bwa.yaml"
	threads: 4
	shell:
		"samtools index -@ {threads} {input}"

rule bigwig:
	input:
		bam = rules.sort.output,
		bai = rules.index.output
	output:
		"data/bigwig/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	threads: 16
	shell:
		"bamCoverage -b {input.bam} -o {output} -p {threads} --skipNAs --binSize 10 --smoothLength 50 --normalizeUsing CPM"
