rule deseq_init:
	input:
		catalog = "data/counts_table.txt"
	output:
		norm_counts = "data/de/norm_counts.txt",
		log_norm_counts = "data/de/log_norm_counts.txt",
		pca = "data/de/sample_PCA.png"
		# heatmap = "data/de/sample_heatmap.png"
	params:
		metadata = config["metadata_file"]
	conda:
		"../envs/deseq2.yaml"
	script:
		"../scripts/deseq2.R"