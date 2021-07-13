rule deseq2:
	input:
		catalog = "data/counts/counts_table.txt",
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
		"../envs/deseq2.yaml"
	threads: 8
	log:
		"data/logs/deseq2.log"
	script:
		"../scripts/deseq2.R"
# normalize by reads in peaks

rule diffbind:
        input:
                consensus_peaks = "data/macs2/consensus_peaks.bed",
                metadata = config["DIFFBIND_CONFIG"]
        params:
                padj_cutoff = config["padj_cutoff"]
        output:
                directory("data/diffbind")
        conda:
                "../envs/diffbind.yaml"
        threads: 8
        script:
                "../scripts/diffbind.R"
# normalize by entire sequencing depth

rule HOMER:
	input:
		rules.deseq2.output.contrast_combinations
	output:
		directory("data/homer")
	params:
		genome = config["GENOME"]
	log:
		"data/logs/homer.log"
	conda:
		"../envs/homer.yaml"
	threads: 8
	shell:
		"bash scripts/homer.sh -i {input} -g {params.genome}"
