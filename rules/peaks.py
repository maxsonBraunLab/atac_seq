rule macs2:
	input:
		bam = "data/shift/{sample}.shifted.filtered.markd.sorted.bam",
		bai = "data/shift/{sample}.shifted.filtered.markd.sorted.bam.bai"
	output:
		"data/macs2/{sample}/{sample}_peaks.broadPeak"
	conda:
		"../envs/macs2.yaml"
	params:
		outdir = "data/macs2/{sample}",
		genome_size = config["GSIZE"]
	log:
		"data/logs/macs2_{sample}.log"
	shell:
		"macs2 callpeak -t {input.bam} -n {wildcards.sample} "
		"--format BAMPE --gsize {params.genome_size} "
		"--outdir {params.outdir} --broad > {log} 2>&1"

rule consensus:
	input:
		expand("data/macs2/{sample}/{sample}_peaks.broadPeak", sample = SAMPLES)
	output:
		"data/macs2/consensus_peaks.bed"
	params:
		n_intersects = config["N_INTERSECTS"]
	conda:
		"../envs/consensus.yaml"
	script:
		"../scripts/consensus_peaks.R"

rule counts:
	input:
		consensus = rules.consensus.output,
		sample = "data/shift/{sample}.shifted.filtered.markd.sorted.bam"
	output:
		"data/multicov/{sample}.txt"
	conda:
		"../envs/bedtools.yaml"
	shell:
		"bedtools multicov -bams {input.sample} -bed {input.consensus} > {output}"

rule counts_table:
	input:
		expand("data/multicov/{sample}.txt", sample = SAMPLES)
	output:
		"data/counts/counts_table.txt"
	run:
		dfs = []
		for file in list(input):
			sample_name = os.path.basename(file).split('.')[0]
			dfs.append( pd.read_csv(file, sep = "\t", names = ["chr", "start", "end", "name", sample_name]) )
		df = pd.concat(dfs, axis = 1)
		df = df.loc[:,~df.columns.duplicated()]
		df.to_csv(str(output), header = True, index = False, sep = "\t")

rule intervene:
	input:
		expand("data/macs2/{sample}/{sample}_peaks.broadPeak", sample = SAMPLES)
	output:
		"data/intervene/intervene_results/amleto_intersect_venn.png"
	conda:
		"../envs/intervene.yaml"
	shell:
		"bash scripts/intervene.sh"