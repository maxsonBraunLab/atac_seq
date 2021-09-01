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
		"bwa mem -t {threads} {config[GENOME]} {input.r1} {input.r2} 2>{log} | samtools sort -@ {threads} > {output}"

rule filter:
	input:
		rules.bwa.output
	output:
		temp("data/filter/{sample}.filtered.sorted.bam")
	conda:
		"../envs/bwa.yaml"
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
		"../envs/sambamba.yaml"
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
		"../envs/bedtools.yaml"
	threads: 4
	shell:
		"bedtools intersect -v -ubam -abam {input.bam} -b {input.banlist} | samtools sort -@ {threads} > {output[0]}; samtools index {output[0]}"

rule bigwig:
	input:
		bam = rules.banlist.output[0],
		bai = rules.banlist.output[1]
	output:
		"data/bigwig/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	threads: 16
	shell:
		"bamCoverage -b {input.bam} -o {output} -p {threads} --skipNAs --binSize 10 --smoothLength 50 --normalizeUsing CPM"
