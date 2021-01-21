rule makeBigwig:
	input:
		"samples/bamfiles/{sample}_rmChrM_dedup_quality_shifted.bam"
	output:
		"data/bigwigs/{sample}_tracks.bw",
	conda:
		"../envs/deeptools.yaml"
	threads: 16
	message:
		"""--- making bigwig---"""
	shell:
		"bamCoverage --bam {input} -o {output} --p {threads} --skipNAs --binSize 10 --smoothLength 50 --normalizeUsing CPM"

#ATAC sequencing involves a transposase that slightly shifts the reads when 
#it is used. In order to correct this we need to 5, 4 offset the reads so that they line up with the actual accessibility locations
rule shift_reads:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam",
		index = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam.bai",
	output:
		tmp = temp("samples/bamfiles/{sample}_tmp.bam"),
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shifted.bam"
	params:
	conda:
		"../envs/deeptools.yaml"
	threads: 4 
	message:
		"""--- shifting reads---"""
	shell:
		"""
		alignmentSieve --bam {input.bamfile} --outFile {output.tmp} -p {threads} --ATACshift
		samtools sort -@ {threads} -o {output.bamfile} {output.tmp}
		samtools index -@ {threads} {output.bamfile}
		"""

rule index:
	input:
		deduplicated = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam"
	output:
		indexed = temp("samples/bamfiles/{sample}_rmChrM_dedup_quality.bam.bai")
	conda:
		"../envs/sambamba.yaml"
	threads: 4 
	message:
		"""--- indexing reads---"""
	shell:
		"sambamba index -t {threads} {input.deduplicated} {output.indexed}"

rule align_stats:
	input:
		expand("samples/sample_stats/{sample}_read_depths.csv", sample = SAMPLE)
	output:
		"samples/sample_stats/align_stats.txt"
	shell:
		"cat {input} > {output}"

#Filter out reads that match the quality flags that we want in samtools
rule hq_mapped_reads:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup.bam"
	output:
		quality = temp("samples/bamfiles/{sample}_rmChrM_dedup_quality.bam"),
	params:
		readsfile = "samples/sample_stats/{sample}_read_depths.csv"
	conda:
		"../envs/samtools.yaml"
	threads: 4
	message:
		"""--- select high-quality mapped reads ---"""
	shell:
		"""
		samtools view -h -@ {threads} -b -F 1804 {input.bamfile} > {output.quality}
		samtools view -@ {threads} -c {output.quality}| xargs -I{{}} echo "{wildcards.sample},satisfy_quality,"{{}}$'\n'  >> {params.readsfile}
		"""
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2

#deduplicate the reads that we have. This will get rid of potential PCR duplicates
#After this step you might have a large difference in reads if some sample have a high percentage of duplicates 
rule deduplicate_reads:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM.bam",
	output:
		temp("samples/bamfiles/{sample}_rmChrM_dedup.bam.bai"),
		deduplicated = temp("samples/bamfiles/{sample}_rmChrM_dedup.bam")
	params:
		readsfile = "samples/sample_stats/{sample}_read_depths.csv",
	conda:
		"../envs/sambamba.yaml"
	threads: 4
	message:
		"""--- deduplicating reads---"""
	shell:
		 """
		 sambamba markdup -r -t {threads} {input.bamfile} {output.deduplicated}
		 samtools view -@ {threads} -c {output.deduplicated}| xargs -I{{}} echo "{wildcards.sample},no_duplicates,"{{}}$'\n'  >> {params.readsfile}
		 """

#Remove the mitochondrial reads This requires that the mitochondrial reads have the  name chrM so that you can grep them with an inverted search
#This works for our current aligner and references, but if you get zero mitochondrial reads filtered out its possible that the mitochondrial chromosome
#has a name that is different from "chrM"
rule rm_mito_reads:
	input:
		bamfile = "samples/align/sorted/{sample}_sorted.bam", #change this to _merged if merging lanes
	output:
		bamfile = temp("samples/bamfiles/{sample}_rmChrM.bam"),
		readsfile = "samples/sample_stats/{sample}_read_depths.csv"
	conda:
		"../envs/samtools.yaml"
	message:
		"""--- Removing mitochondrial reads---"""
	threads: 4 
	shell:
		"""
		samtools view -@ {threads} -c {input.bamfile}| xargs -I{{}} echo "{wildcards.sample},starting_reads,"{{}}$'\n' >> {output.readsfile}
		samtools view -h -@ {threads} {input.bamfile}| grep -v chrM| samtools sort -@ {threads} -O BAM > {output.bamfile}
		samtools view -@ {threads} -c {output.bamfile}| xargs -I{{}} echo "{wildcards.sample},no_ChrM,"{{}}$'\n' >> {output.readsfile}
		"""

