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

#This is the location for our merged narrow peaks file
np_file <- args[1]
blacklist_file <- args[2]
genome_name <- "hg38"
#genome_name <- args[5]

presence_in_samples <- as.numeric(args[3])
#set the name that we want for the final out file
bedfile <- args[4]



#This command reads in the merged broadnarrowpeak files and adds the information about the genome
#we used 
#peaks <- read_narrowPeaks(np_file, genome_info = genome_name)
peaks <- read.delim(np_file, header=FALSE)
colnames(peaks) <- c("seqnames", "start", "end", "name", "score", "no_strand", "signalValue", "pValue","qvalue")
peaks <- as_granges(peaks) %>% set_genome_info(genome = genome_name)

#now lets read in the blacklist file so we can remove blacklisted regions
hg38.blacklist.v2 <- read.delim(blacklist_file, header=FALSE)
colnames(hg38.blacklist.v2) <- c("seqnames","start","end","reason")
blacklist <- as_granges(hg38.blacklist.v2) %>% set_genome_info(genome = genome_name)

#This command computes the coverage at each location on the genome
#eg. how many overlapping peaks exist at this particular location
cvg <- peaks %>% 
  compute_coverage()

#To create our peak catalog we filter so that every location is overlapped by at least two peaks
#We then collate together peaks that are right next to each other
peaks_catalog <- cvg %>% 
  filter(score >= presence_in_samples)


#filter out only the chromosomes that we want and that we care about
peaks_catalog <- peaks_catalog %>% filter(., str_detect(seqnames,"^chr[\\dXY]+$"))

#now we filter out everything that is overlapping with the blacklist file
peaks_catalog <- peaks_catalog %>% filter_by_non_overlaps(blacklist)

#stretch the ranges by one so that the reduced ranges will be added together
peaks_catalog_stretched <- peaks_catalog %>%         
  anchor_start() %>%
  stretch(1)

#reduce the ranges to merge them
peaks_merged <- peaks_catalog %>%
  reduce_ranges()
# Then to get rid of the tiny amount at the end
# ranges %>% anchor_start() %>% stretch(-1)


write_bed(peaks_merged, bedfile)

