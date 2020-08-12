args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
samples <- read.csv(args[1], header=FALSE)

#V1","starting_reads","satisfy_quality","reads_in_peaks","percent_mapped","number_of_mitochondrial","number_of_duplicates","frip"
dataframe <- pivot_wider(samples, names_from=V2, values_from = V3)
out_df <- data.frame(matrix(ncol = 1 , nrow = nrow(dataframe)))
out_df$satisfy_quality <- dataframe$satisfy_quality
out_df$starting_reads <- dataframe$starting_reads
out_df$number_of_peaks <- dataframe$unfiltered_peaks
out_df$number_of_mito <- dataframe$starting_reads - dataframe$no_ChrM
out_df$number_of_duplicates <- dataframe$no_ChrM - dataframe$no_duplicates
out_df$reads_in_peaks <- dataframe$reads_in_peaks
out_df$percent_mapped <- dataframe$mapping_reads / dataframe$starting_reads
out_df$frip <- dataframe$reads_in_peaks / dataframe$satisfy_quality
out_df$sample_name <- dataframe$V1
out_df$matrix.ncol...1..nrow...nrow.dataframe..<- NULL
write.table(out_df, row.names=F, file=args[2], sep='\t', quote=F)


