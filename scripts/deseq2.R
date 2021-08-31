sink(file("data/logs/deseq2.log", open = "wt"), type = "message")

library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(tibble)

parallel <- FALSE
if (snakemake@threads > 1) {
	library("BiocParallel")
	# setup parallelization
	register(MulticoreParam(snakemake@threads))
	parallel <- TRUE
}

# import snakemake variables --------------------------------------------------
catalog <- snakemake@input[['catalog']]
md <- snakemake@input[['metadata']]

print(paste("Peaks catalog: ", catalog))
print(paste("Metadata file: ", md))

# import data ---------------------------------------------------------
catalog <- read.table(catalog, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
md <- read.table(md, sep = "\t", header = T, stringsAsFactors = F)
consensus_intervals <- catalog[,1:4]

print(head(catalog))
print(head(md))

counts <- catalog[,md$SampleID]
rownames(counts) <- catalog$name
rownames(md) <- md$SampleID

# deseq2 ----------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = md, design = ~Condition)
dds <- DESeq(dds, parallel = parallel)
res <- results(dds)

# export normalized data
NormCounts <-counts(dds, normalized=TRUE)
NormCounts <- 0.00001+(NormCounts) #this makes it so there are no 0 values#
LogNormCounts <- log2(NormCounts)
LogNormCounts <- as.data.frame(LogNormCounts)
LogNormCounts$ID <- rownames(LogNormCounts)

write.csv(NormCounts, snakemake@output[['norm_counts']], quote = F)
write.csv(LogNormCounts, snakemake@output[['log_norm_counts']], quote = F)
print("Exported normalized counts tables")

# clustering - PCA
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "Condition") + geom_label_repel(aes(label=name))
ggsave(snakemake@output[['pca']], pca, width = 16, height = 9, dpi = 300, units = "in")

# define all contrasts --------------------------------------------------------
contrast_combinations <- noquote(t(combn(unique(md$Condition), 2)))
colnames(contrast_combinations) <- c("c1", "c2")
print("Combinations of contrasts to assess")
print(contrast_combinations)
write.table(contrast_combinations, snakemake@output[["contrast_combinations"]], sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

de_stats <- data.frame()
all_sig_intervals <- data.frame()

if (!dir.exists("data/deseq2")) {
	dir.create("data/deseq2")
}

# DE results per contrast combination -----------------------------------------

for (i in 1:nrow(contrast_combinations)) {

	# define contrasts
	c1 = contrast_combinations[i, "c1"]
	c2 = contrast_combinations[i, "c2"]

	# define outputs + make parent dir.
	contrast_name <- paste0(c1, "-vs-", c2)
	output_prefix <- paste0("data/deseq2/", contrast_name, "/")
	if (!dir.exists(output_prefix)) {dir.create(output_prefix)}

	all_file <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-all.txt") # contain all peaks without NA padj
	sig_file <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-sig.txt") # all peaks padj < padj_cutoff
	up_file <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-up_sig.txt") # all peaks with logFC > 0
	down_file <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-down_sig.txt") # all peaks with logFC < 0
	up_intervals <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-up_sig.bed") # all peaks with logFC > 0 in bed format
	down_intervals <- paste0("data/deseq2/", contrast_name, "/", contrast_name, "-down_sig.bed") # all peaks with logFC < 0 in bed format

	# extract contrasts from DESeq object, identify up / down results and intervals.
	temp_res <- results(dds, contrast = c("Condition", c1, c2), cooksCutoff = FALSE)
	temp_res_df <- temp_res %>% 
			as.data.frame() %>% 
			rownames_to_column("name") %>% 
			arrange(padj) # contains all peaks and DE info.

	all_sig <- temp_res_df %>% filter(padj <= snakemake@params[['padj_cutoff']])
	upregulated_sig <- temp_res_df %>% filter(log2FoldChange > 0 & padj <= snakemake@params[['padj_cutoff']])
	downregulated_sig <- temp_res_df %>% filter(log2FoldChange < 0 & padj <= snakemake@params[['padj_cutoff']])
	
	upregulated_intervals <- upregulated_sig %>% 
		inner_join(x = ., y = consensus_intervals, by = "name") %>% 
		select(chr, start, end, name)

	downregulated_intervals <- downregulated_sig %>%
		inner_join(x = ., y = consensus_intervals, by = "name") %>% 
		select(chr, start, end, name)

	temp_res_df <- temp_res_df %>% rename(name = "name") # intervals should be sorted by padj
	temp_sig_intervals <- rbind(upregulated_sig, downregulated_sig)
	all_sig_intervals <- rbind(all_sig_intervals, temp_sig_intervals)

	write.table(temp_res_df, all_file, quote = FALSE, sep = "\t", row.names = F)
	write.table(all_sig, sig_file, quote = FALSE, sep = "\t", row.names = F)
	write.table(upregulated_sig, up_file, quote = FALSE, sep = "\t", row.names = F)
	write.table(downregulated_sig, down_file, quote = FALSE, sep = "\t", row.names = F)
	write.table(upregulated_intervals, up_intervals, quote = FALSE, sep = "\t", row.names = F, col.names = FALSE)
	write.table(downregulated_intervals, down_intervals, quote = FALSE, sep = "\t", row.names = F, col.names = FALSE)

	temp_stats <- data.frame(cond1 = c1, cond2 = c2, tot_peaks = nrow(temp_res_df), sig_peaks = nrow(all_sig), sig_up = nrow(upregulated_sig), sig_down = nrow(downregulated_sig))
	de_stats <- rbind(de_stats, temp_stats)

	print(paste("Results for", contrast_name, "has been exported"))
	print(paste(c1, c2))
	print("Significant results")
	print(all_sig)
}

all_sig_intervals <- all_sig_intervals %>%
	inner_join(x = ., y = consensus_intervals, by = "name") %>%
	select(chr, start, end, name) %>%
	rename(name = "name")

write.table(de_stats, snakemake@output[['stats']], quote = FALSE, sep = "\t", row.names = F, col.names = T)
write.table(all_sig_intervals, snakemake@output[['all_sig_intervals']], quote = FALSE, sep = "\t", row.names = F, col.names = F)