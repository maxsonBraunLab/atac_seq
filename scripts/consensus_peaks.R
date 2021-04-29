library(dplyr)
library(stringr)
library(GenomicRanges)
library(plyranges)
library(tidyr)

# for each condition, read in peaks containing all samples.
list_of_peak_files <- lapply(unlist(snakemake@input),
					read.delim,
					header = FALSE,
					col.names = c("seqnames", "start",
						"end", "name",
						"score", "strand",
						"signalValue", "pval",
						"qval", "peak"))
all_peaks_df <- bind_rows(list_of_peak_files)

# define condition and replicate.
all_peaks_df <- all_peaks_df %>%
	separate(name, "condrep", sep = "_", remove=FALSE, extra = "drop") %>%
	mutate(rep = str_sub(condrep, start = -1)) %>%
	mutate(cond = str_sub(condrep, start = 1, end = -2)) %>%
	mutate(condrep = NULL) %>%
	mutate(strand = "*")

# split the df into list of GRanges by condition. compute coverage for each condition.
peaks_by_condition <- makeGRangesListFromDataFrame(all_peaks_df, split.field = "cond")

# find peaks present in this many or more samples in one condition.
presence_in_samples <- snakemake@params$n_intersects

consensus_peaks <- GRangesList()
for (condition in names(peaks_by_condition)) {
	overlaps_by_condition <- compute_coverage(peaks_by_condition[condition]) %>%
		filter(score >= presence_in_samples) %>%
		reduce_ranges()

	if (length(overlaps_by_condition) == 0) {
		next
	} else {
		seqlengths(overlaps_by_condition) <- rep(NA,length(seqlengths(overlaps_by_condition))) # set seqlengths as NA
		consensus_peaks[[condition]] <- overlaps_by_condition
	}
}

# You can see which conditions contributed to the consensus peaks with the following!
# unlist(consensus_peaks)

# merge overlapping peaks across conds, turn list of GRanges into df, rm blacklist regions, export.
consensus_df <- unlist(consensus_peaks) %>% reduce_ranges()
consensus_df <- as.data.frame(consensus_df, row.names = NULL)
consensus_df <- consensus_df %>% 
					select(seqnames, start, end) %>%
					mutate(name = paste("consensus_peak", row_number(), sep = "_")) %>%
					filter(start < end)

write.table(x = consensus_df,
	file = snakemake@output[[1]], sep = "\t",
	row.names = FALSE, col.names = FALSE,
	quote = FALSE)