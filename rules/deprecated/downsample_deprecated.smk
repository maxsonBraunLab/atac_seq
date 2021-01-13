#This command creates the bigwig files to view from after the reads have been fully preprocessed
#This command can take a while and eventually needs to be replaced with the smk process in the create_merged_bigwig.smk file
rule makeBigwig:
	input:
		quality = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam",
	output:
		bigwig = "data/bigwigs/{sample}_tracks_5window_smooth.bw",
		bigwig_rough = "data/bigwigs/{sample}_tracks_5window_rough.bw",
	conda:
		"../envs/deeptools.yaml"
	threads: 16
	message:
		"""--- making bigwig---"""
	shell:
		"""
		bamCoverage --bam {input.quality} -o {output.bigwig_rough} --numberOfProcessors {threads} --skipNAs --binSize 5 
		bamCoverage --bam {input.quality} -o {output.bigwig} --numberOfProcessors {threads} --skipNAs --binSize 5 --smoothLength 15  
		"""

#This rule downsapmles the bam files by taking an amount to downsample a file from the file downsample.list
#Then this command uses the tool samtools -s to downsample the file to the correct percentage
#currently this is not exact down to the number of reads but give me bamfiles within around 1% of each other
#This rule uses a lot of complicated bash trickery and should be made more clear In line 4 I first grep the file to check to make sure that
#the sample is in the downsample list, there is only one sample that is not and this is hte sample that will not be downsample that we are downsampling everything to
#So if that is the sample we simply copy it to the _downsampled.bam, otherwise we grep out the sample name and take the second column where the sample percentage is and downsample the sample that amount
rule downSample:
	input:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads.bam",
		downsample_list = "samples/bamfiles/downsample_list.txt",
	output:
		bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam",
		index = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam.bai",
		input_index = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam.bai",
		input_sort = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam",
	conda:
		"../envs/sambamba.yaml"
	threads: 4 
	message:
		"""---normalizing reads---"""
	shell:
		"""
		sambamba sort -t {threads} {input.bamfile} -o {output.input_sort}
		sambamba index -t {threads} {output.input_sort} {output.input_index}
		export PATH=$PATH:/opt/installed/samtools-1.6/bin/
		if [[ $(grep {wildcards.sample} {input.downsample_list}) ]]; then samtools view -@ 4 -b -h -s 1.$(grep {wildcards.sample} {input.downsample_list}| awk '{{print $2}}') {output.input_sort} -o {output.bamfile}; else cp {output.input_sort} {output.bamfile}; fi 
		sambamba index -t {threads} {output.bamfile} {output.index}
		"""

rule createDownsample:
	input:
		bamfile = expand("samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads.bam", sample=SAMPLE),
		readsfile = expand("samples/sample_stats/{sample}_read_depths.csv", sample=SAMPLE),
	output:
		downsample_list = "samples/bamfiles/downsample_list.txt",
		filtered_list = "samples/bamfiles/sample_downsample_information.txt",
	message:
		"""--calculating normalization factor to subsample by---"""
	run:
		import pandas as pd
		dataframes = []
		for file in input.readsfile:
			dataframes.append(pd.read_csv(file, names=["name","description","reads"]))
		merged_dataframe = pd.concat(dataframes, ignore_index=True)
		print(merged_dataframe)        
		filtered = merged_dataframe[merged_dataframe.description=="satisfy_quality"] #ensure that we are normalizing to the final processed reads
		print(filtered)
		filtered = filtered.drop_duplicates() #its possible that some values have been written more than once
		print(filtered)                       #we want to drop these values as they are identical rows and are not needed
		filtered['downsample'] = pd.to_numeric(filtered['reads'].min()) / pd.to_numeric(filtered['reads']) #what percentage of the smallest sample to downsample
		print(filtered)  # each of the sample to then print the result
		filtered.to_csv(output.filtered_list, header=True, index=False, sep='\t')
		out_data = filtered.loc[:,["name", "downsample"]]
		normalizing_sample = out_data[out_data.downsample == 1] # which one is hte normalizing sample
		print(normalizing_sample)
		out_data = out_data[out_data.downsample < 1] #only keep sample that need to be downsampled and skip the ones that wont be
		out_data["downsample"] = (round(out_data["downsample"], 8)*100000000).astype(int) #round the sample to the nearest number, I will add the . later
		out_data.to_csv(output.downsample_list, header=False, index=False, sep='\t') #write everything out to a csv


