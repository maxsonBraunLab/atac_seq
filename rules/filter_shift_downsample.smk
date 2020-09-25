#ATAC sequencing involves a transposase that slightly shifts the reads when 
#it is used. In order to correct this we need to 5, 4 offset the reads so that they line up with the actual accessibility locations
rule shiftReads:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam",
		index = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam.bai",
	output:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads.bam",
	conda:
		"../envs/deeptools.yaml"
	threads: 4 
	message:
		"""--- shifting reads---"""
	shell:
		"""
		alignmentSieve --bam {input.bamfile} --outFile {output.bamfile} --numberOfProcessors {threads} --ATACshift
		"""

#many commands require that there is index for the bamfile that you are using. TODO: this command should
#eventually be pulled together into the command that needs this index which is hte shifting rule  above
rule index:
	input:
		deduplicated = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam",
		full_bam = "samples/align/sorted/{sample}_sorted.bam",
	output:
		indexed = "samples/bamfiles/{sample}_rmChrM_dedup_quality.bam.bai",
		indexed_full = "samples/bamfiles/{sample}_sorted.bam.bai",
	conda:
		"../envs/sambamba.yaml"
	threads: 4 
	message:
		"""--- indexing reads---"""
	shell:
		"""
		sambamba index -t {threads} {input.full_bam} {output.indexed_full}
		sambamba index -t {threads} {input.deduplicated} {output.indexed}
		
		"""

#Filter out reads that match the quality flags that we want in samtools
rule samtoolsQuality:
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
		"""--- deduplicating reads---"""
	shell:
		"""
		samtools view -h -@ 4 -b -F 1804 {input.bamfile} > {output.quality}
		samtools view -@ 4 -c {output.quality}| xargs -I{{}} echo "{wildcards.sample},satisfy_quality,"{{}}$'\n'  >> {params.readsfile}
		"""
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2

#deduplicate the reads that we have. This will get rid of potential PCR duplicates
#After this step you might have a large difference in reads if some sample have a high percentage of duplicates 
rule dedup:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM.bam",
	output:
		deduplicated = temp("samples/bamfiles/{sample}_rmChrM_dedup.bam"),
	params:
		readsfile = "samples/sample_stats/{sample}_read_depths.csv",
	conda:
		"../envs/sambamba.yaml"
	threads: 4
	message:
		"""--- deduplicating reads---"""
	shell:
		 """
		 sambamba markdup -t {threads} {input.bamfile} {output.deduplicated}
		 samtools view -@ 4 -F 1024 -c {output.deduplicated}| xargs -I{{}} echo "{wildcards.sample},no_duplicates,"{{}}$'\n'  >> {params.readsfile}
		 """

#Remove the mitochondrial reads This requires that the mitochondrial reads have the  name chrM so that you can grep them with an inverted search
#This works for our current aligner and references, but if you get zero mitochondrial reads filtered out its possible that the mitochondrial chromosome
#has a name that is different from "chrM"
rule removeMitochondrial:
	input:
		bamfile = "samples/align/sorted/{sample}_sorted.bam", #change this to _merged if merging lanes
	output:
		bamfile = temp("samples/bamfiles/{sample}_rmChrM.bam"),
		readsfile = "samples/sample_stats/{sample}_read_depths.csv"
	conda:
		"../envs/samtools.yaml"
	message:
		"""--- Removing mitochondrial reads---"""
	log:
		totalreads = "samples/bamfiles/{sample}_readswmito.txt"
	threads: 4 
	shell:
		"""
		samtools view -@ {threads} -c {input.bamfile}| xargs -I{{}} echo "{wildcards.sample},starting_reads,"{{}}$'\n' >> {output.readsfile}
		samtools view -h -@ {threads} {input.bamfile}| grep -v chrM| samtools sort -@ {threads} -O BAM > {output.bamfile}
		samtools view -@ {threads} -c {output.bamfile}| xargs -I{{}} echo "{wildcards.sample},no_ChrM,"{{}}$'\n' >> {output.readsfile}
		"""
