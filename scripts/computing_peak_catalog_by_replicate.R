#Script will first create a coverage graph then filter out antying that does not have coverage of at least 3
#After that will restitch together using reduce ranges
#Finally we will count the number of overlaps that we are getting from bam files over the ranges that we have decided to use

#we used the hg38 genome vs the hg19 genome because it was more up to date and included mitochondrial
# allowing us to quantify how many actual mitochondrial reads we got that were mapped
#first we are going to get the arguments that are passed in from the command line
args = commandArgs(trailingOnly=TRUE)
#The format for the input will be, Rscript this_script_name narrowpeaksfile blacklistfile presence_in_samples_no outbedfile


library(dplyr)
library(plyranges)
library(stringr)
library(purrr)
library(tidyr)
library(Rsamtools)
library(parallel)

#This is the location for our merged narrow peaks file
np_file <- args[1]
blacklist_file <- args[2]

#genome_name <- "hg38"
genome_name <- args[5]

#presence_in_samples <- 2
presence_in_samples <- as.numeric(args[3])

#set the name that we want for the final out file
bedfile <- args[4]

#read metadata file
metadata_file <- args[6]
replace_expression <- args[7]

options(srapply_fapply="parallel", mc.cores = 8)

#This command reads in the merged broadnarrowpeak files and adds the information about the genome
#we used 
#peaks <- read_narrowPeaks(np_file, genome_info = genome_name)
peaks <- read.delim(np_file, header=FALSE)

#lets load the peaks file and convert it to a granges object
colnames(peaks) <- c("seqnames", "start", "end", "name", "score", "no_strand", "signalValue", "pValue","qvalue")
peaks <- as_granges(peaks) %>% set_genome_info(genome = genome_name)

#now lets load the metadata file
metadata <- read.delim(metadata_file, header=TRUE)
print("metadata imported")
#now lets read in the blacklist file so we can remove blacklisted regions
hg38.blacklist.v2 <- read.delim(blacklist_file, header=FALSE)
colnames(hg38.blacklist.v2) <- c("seqnames","start","end","reason")
blacklist <- as_granges(hg38.blacklist.v2) %>% set_genome_info(genome = genome_name)

#filter out only the chromosomes that we want and that we care about
#This filtering command keeps all chromosomes that are numbered with nothing after the
#number, as well as the X and Y chromosomes
filtered_peaks <- peaks %>% filter(., str_detect(seqnames,"^chr[\\dXY]+$"))

# create extra col to indicate treatment + define replicate
filtered_peaks <- as.data.frame(filtered_peaks) %>% separate(name, into = c("treatment"), sep = "_", remove=FALSE)
rep_df <- as.data.frame(str_split_fixed(filtered_peaks$name, "_", 3))
filtered_peaks$replicate <- paste(rep_df$V1, rep_df$V2, sep = "_")
filtered_peaks <- as_granges(filtered_peaks)

#now we are going to make a list of all the treatments by group and essentially apply reduce
#reduce doesn't work here so instead we are using a for loop
# for each condition, find consensus peaks.
peaks_list <- filtered_peaks %>% split(., filtered_peaks$treatment)
catalog <- NULL
for (item in 1:length(peaks_list)) {

    coverage <- compute_coverage(peaks_list[item]) %>%
       filter(score >= presence_in_samples) %>%
       reduce_ranges() %>%
       filter_by_non_overlaps(blacklist)

    if (length(catalog)==0) {
       catalog <- coverage 
    } else {
       seqlengths(catalog) <- rep(NA,length(seqlengths(catalog))) #this fixes an error where it recalculates seq lengths giving different lengths
       seqlengths(coverage) <- seqlengths(catalog)                # and then won't let you create a union between them
       catalog <- GenomicRanges::union(catalog, coverage)
    }
}

catalog <- catalog %>% as.data.frame() %>% mutate(name = paste0("consensus_peak_",row_number())) %>% as_granges()
write_bed(catalog, bedfile)

# for each replicate, write stats ----------------------------------------------------------------------------
rep_split <- filtered_peaks %>% split(., filtered_peaks$replicate)
stats_list <- lapply(rep_split, function(x) {
    replicate <- toString(unique(x$replicate))
    treatment <- unique(x$treatment)
    peaks_in_replicate <- as.integer(nrow(as.data.frame(x)))
    peaks_in_consensus <- nrow(as.data.frame(  join_overlap_inner(x, catalog)  ))
    number_consensus_peaks <- nrow(as.data.frame(catalog))
    temp_stats <- data.frame(
        treatment=treatment,
        rep = replicate,
        tot_peaks = peaks_in_replicate,
        peaks_in_consensus=peaks_in_consensus,
        tot_consensus_peaks=number_consensus_peaks)
})
consensus_stats <- do.call(rbind, stats_list)
save.image()

# lastly, add FRCC metric = fraction of reads in consensus catalog
catalog_df <- as.data.frame(catalog)
for (i in 1:nrow(consensus_stats)) {
  # define input files
  sample <- consensus_stats[i, "rep"]
  in_file <- paste0("samples/bamfiles/filtered/", sample, "_rmChrM_dedup_quality_shiftedReads_sorted.bam")
  bamFile <- BamFile(in_file)
  # in regions is catalog variable

  # count reads in consensus intervals
  gr <- GRanges( seqnames = catalog_df$seqnames, ranges = IRanges(start = c(catalog_df$start), end = c(catalog_df$end) ))
  counted <- countBam(bamFile, param = ScanBamParam(which=gr, what = scanBamWhat() ))
  reads_in_cc <- sum(counted$records)
  total_reads <- countBam(bamFile, param = ScanBamParam(what = "qname"))$records
  frcc <- reads_in_cc / total_reads
  consensus_stats[i, "tot_reads"] <- total_reads
  consensus_stats[i, "reads_in_CC"] <- reads_in_cc
  consensus_stats[i, "%_reads_in_CC"] <- frcc
}

consensus_stats
write.table(consensus_stats, "samples/macs/consensus_stats.txt", quote= FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)