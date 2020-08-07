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


#This command reads in the merged broadnarrowpeak files and adds the information about the genome
#we used 
#peaks <- read_narrowPeaks(np_file, genome_info = genome_name)
peaks <- read.delim(np_file, header=FALSE)

#lets load the peaks file and convert it to a granges object
colnames(peaks) <- c("seqnames", "start", "end", "name", "score", "no_strand", "signalValue", "pValue","qvalue")
peaks <- as_granges(peaks) %>% set_genome_info(genome = genome_name)

#now lets load the metadata file
metadata <- read.delim(metadata_file, header=TRUE)

#now lets read in the blacklist file so we can remove blacklisted regions
hg38.blacklist.v2 <- read.delim(blacklist_file, header=FALSE)
colnames(hg38.blacklist.v2) <- c("seqnames","start","end","reason")
blacklist <- as_granges(hg38.blacklist.v2) %>% set_genome_info(genome = genome_name)

#filter out only the chromosomes that we want and that we care about
#This filtering command keeps all chromosomes that are numbered with nothing after the
#number, as well as the X and Y chromosomes
filtered_peaks <- peaks %>% filter(., str_detect(seqnames,"^chr[\\dXY]+$"))

#now we are going to edit and fix the names on the all_broad_peaks table so we can match them with the sample ID
#replace_expression <- "-ATAC/.*"
filtered_peaks$name <- as.character(filtered_peaks$name) %>% str_replace(., replace_expression,"")


metadata$Condition <- as.character(metadata$Condition)
#match the name of the condition to the name of the cleaned sample ID
filtered_peaks$treatment <- as.character(map(filtered_peaks$name, (function (x) metadata[as.character(metadata$SampleID) == as.character(x),]$Condition)))
print(unique(filtered_peaks$treatment))
print(unique(filtered_peaks$name))

print(unique(metadata$SampleID))
print(unique(metadata$Condition))

#now we are going to make a list of all the treatments by group and essentially apply reduce
#reduce doesn't work here so instead we are using a for loop
peaks_list <- filtered_peaks %>% split(., filtered_peaks$treatment)
print("now for the names of the peaks list")
print(names(peaks_list))

catalog <- NULL
names(peaks_list)
for (item in 1:length(peaks_list)) {

    coverage <- compute_coverage(peaks_list[item]) %>%
       filter(score >= presence_in_samples) %>%
       reduce_ranges() %>%
       filter_by_non_overlaps(blacklist)

    print(head(peaks_list[item]))
    #Now lets get the name of what the condition is
    #since we split into lists by name, the name of the list we are using
    #will be associated with the condition that it is
    treatment_text <- names(peaks_list)[item]
    print(treatment_text)
    write_bed(coverage, paste0(bedfile, "_", treatment_text, ".bed") ) 

    if (length(catalog)==0) {
       catalog <- coverage 
    }else{
       seqlengths(catalog) <- rep(NA,length(seqlengths(catalog))) #this fixes an error where it recalculates seq lengths giving different lengths
       seqlengths(coverage) <- seqlengths(catalog)                # and then won't let you create a union between them
       catalog <- GenomicRanges::union(catalog, coverage)

    }
        
}

cleaned_catalog <- catalog   #%>% filter_by_non_overlaps(blacklist)

named_catalog <- cleaned_catalog %>% as.data.frame() %>% mutate(name = paste0("consensus_peak_",row_number())) %>% as_granges()
named_catalog

write_bed(named_catalog, bedfile)

