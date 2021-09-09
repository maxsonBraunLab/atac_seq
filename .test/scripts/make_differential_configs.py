#!/usr/bin/env python3

import os
import glob
import sys
import pandas as pd

def main():

	deseq2_outfile = "data/config/deseq2_config.tsv"
	diffbind_outfile = "data/config/diffbind_config.tsv"

	samples = sorted(glob.glob("data/raw/*.fastq.gz"))

	if len(samples) == 0:
		print("ERROR: No samples (ending in .fastq.gz format) found in data/raw.")
		sys.exit()

	df = pd.DataFrame(samples, columns = ["file"])
	df["sample"] = df.file.apply(lambda x: os.path.basename(x).split(".")[0].replace("_R1", "").replace("_R2", ""))
	df["condition"] = df["sample"].apply(lambda x: x.split("_")[0])
	df = df[["sample", "condition"]].drop_duplicates()

	# deseq2 output
	df.to_csv('data/config/deseq2_config.tsv', sep = "\t", header = ["SampleID", "Condition"], index = False)

	# diffbind output
	diffbind = df
	diffbind["replicate"] = diffbind["sample"].apply(lambda x: x.split("_")[1])
	diffbind["bam"] = diffbind["sample"].apply(lambda x: "data/banlist/" + x + ".banlist.filtered.rmdup.sorted.bam")
	diffbind["peaks"] = ["data/macs2/consensus_peaks.bed"] * diffbind.shape[0]
	diffbind["peakcaller"] = ["macs2"] * diffbind.shape[0]
	diffbind["peakformat"] = ["bed"] * diffbind.shape[0]
	diffbind.to_csv('data/config/diffbind_config.csv', sep = ",", header = ["SampleID", "Condition", "Replicate", "bamReads", "Peaks", "PeakCaller", "PeakFormat"], index = False)

if __name__ == "__main__": 
	main()