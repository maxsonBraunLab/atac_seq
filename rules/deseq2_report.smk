# automatically assesses all contrast combinations, exports deseq2-normalized counts tables, and PCA plots for all samples
rule deseq2:
	input:
		catalog = "data/counts_table.txt",
		metadata = config["metadata_file"]
	output:
		directory("data/de"),
		norm_counts = "data/de/norm_counts.txt",
		log_norm_counts = "data/de/log_norm_counts.txt",
		pca = "data/de/sample_PCA.png",
		stats = "data/de/de_stats.txt",
		all_sig_intervals = "data/de/all_sig_intervals.bed"
	params:
		significance = config["significance"]
	conda:
		"../envs/deseq2.yaml"
	threads: 8
	script:
		"../scripts/deseq2.R"

rule essential_report:
	input:
		# pipeline QC metrics
		align_stats = "samples/sample_stats/align_stats.txt",
		consensus_stats = "samples/macs/consensus_stats.txt",
		de_stats = "data/de/de_stats.txt"
	output:
		"data/essential_report.html"
	conda:
		"../envs/report.yaml"
	params:
		# user info
		title = config["title"],
		authors = config["authors"],
		# written content
		intro = config["intro"],
		analysis = config["analysis"],
		takeaways = config["takeaways"],
		notes = config["notes"],
		# config files
		contrasts = config["metadata_file"],
		# misc
		frip_folder = directory("samples/frip")
	script:
		"../scripts/essential_report.Rmd"