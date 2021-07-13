library(DiffBind)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(stringr)

if (dir.exists("data/diffbind")) {
	unlink("data/diffbind", recursive = TRUE)
} else {
	dir.create("data/diffbind")
}

consensus_peaks <- GRanges(read.table(snakemake@input[["consensus_peaks"]], col.names=c("seqnames", "start", "end", "name")))

ATAC <-dba(sampleSheet=snakemake@input[["metadata"]])
ATAC <- dba.count(ATAC, peaks=consensus_peaks)
ATAC <- dba.normalize(ATAC, normalize=DBA_NORM_LIB) # very important
ATAC <- dba.contrast(ATAC, minMembers = 2, categories=DBA_CONDITION)
ATAC <- dba.analyze(ATAC)

# loop through all contrasts and export results
dba_meta <- dba.show(ATAC, bContrasts=TRUE)

# save.image()
for (i in 1:nrow(dba_meta)) {

	de_counts = dba_meta[i, "DB.DESeq2"]
	group1 = dba_meta[i, "Group"]
	group2 = dba_meta[i, "Group2"]
	contrast = paste0(group1, "-vs-", group2)

	if (de_counts == 0) {

		print(paste("WARNING:", contrast, "did not have any differential peaks."))
		next

	} else {

		dba_intervals = dba.report(ATAC, contrast = i, th = 0.05) %>%
			as.data.frame

		# rename Conc_ prefix in columns to "log2MeanReads"
		colnames(dba_intervals) <- str_replace(colnames(dba_intervals), "Conc", "log2MeanReads")

		upregulated_intervals <- dba_intervals %>% filter(Fold > 0)
		downregulated_intervals <- dba_intervals %>% filter(Fold < 0)

		# define output files
		all_sig <- paste0("data/diffbind/", contrast, "/", contrast, "-all_sig.txt")
		up_output <- paste0("data/diffbind/", contrast, "/", contrast, "-up-", "0.05.txt")
		dn_output <- paste0("data/diffbind/", contrast, "/", contrast, "-dn-", "0.05.txt")
		up_bed <- paste0("data/diffbind/", contrast, "/", contrast, "-up-", "0.05.bed")
		dn_bed <- paste0("data/diffbind/", contrast, "/", contrast, "-dn-", "0.05.bed")
		dir.create(paste0("data/diffbind/", contrast), showWarnings = FALSE, recursive = TRUE)

		dba_intervals %>% write.table(., all_sig, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		upregulated_intervals %>% write.table(., up_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		downregulated_intervals %>% write.table(., dn_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		upregulated_intervals %>% select(seqnames, start, end) %>% write.table(., up_bed, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		upregulated_intervals %>% select(seqnames, start, end) %>% write.table(., dn_bed, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
}
