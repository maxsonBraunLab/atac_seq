
#This command creates the bigwig files to view from after the reads have been fully preprocessed
#This command can take a while and eventually needs to be replaced with the smk process in the create_merged_bigwig.smk file
rule makeBigwig:
    input:
        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam",
    output:
        bigwig = sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_smooth.bw",
        bigwig_rough = sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_rough.bw",

    conda:
        "../envs/deeptools.yaml"
    threads: 8 
    message:
        """--- making bigwig---"""
    shell:
        """
        bamCoverage --bam {input.quality} -o {output.bigwig_rough} --numberOfProcessors {threads} --skipNAs --binSize 5 
        bamCoverage --bam {input.quality} -o {output.bigwig} --numberOfProcessors {threads} --skipNAs --binSize 5 --smoothLength 15  
        """
#rule bamcoverage bedtools multicov [OPTIONS] -bams BAM1 BAM2 BAM3 ... BAMn -bed  <BED/GFF/VCF>


#This rule downsapmles the bam files by taking an amount to downsample a file from the file downsample.list
#Then this command uses the tool samtools -s to downsample the file to the correct percentage
#currently this is not exact down to the number of reads but give me bamfiles within around 1% of each other
#This rule uses a lot of complicated bash trickery and should be made more clear In line 4 I first grep the file to check to make sure that
#the sample is in the downsample list, there is only one sample that is not and this is hte sample that will not be downsample that we are downsampling everything to
#So if that is the sample we simply copy it to the _downsampled.bam, otherwise we grep out the sample name and take the second column where the sample percentage is and downsample the sample that amount
rule downSample:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
        downsample_list = sample_work_path + "/bamfiles/downsample_list.txt",
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam",
        index = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam.bai",
        input_index = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam.bai",
        input_sort = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam",
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
        #if [[ $(grep {wildcards.merged_sample} {input.downsample_list}) ]]; then sambamba view --subsample=0.$(grep {wildcards.merged_sample} {input.downsample_list}| awk '{{print $2}}') --subsampling-seed=2 {output.input_sort} -o {output.bamfile}; else cp {output.input_sort} {output.bamfile}; fi
        if [[ $(grep {wildcards.merged_sample} {input.downsample_list}) ]]; then samtools view -@ 4 -b -h -s 1.$(grep {wildcards.merged_sample} {input.downsample_list}| awk '{{print $2}}') {output.input_sort} -o {output.bamfile}; else cp {output.input_sort} {output.bamfile}; fi 
        sambamba index -t {threads} {output.bamfile} {output.index}
        """

rule createDownsample:
    input:
        bamfile = expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam", merged_sample=MERGED_SAMPLES),
        readsfile = expand(sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv", merged_sample=MERGED_SAMPLES),
    output:
        downsample_list = sample_work_path + "/bamfiles/downsample_list.txt",
        filtered_list = sample_work_path + "/bamfiles/samples_downsample_information.txt",
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
        print(filtered)  # each of the samples to then print the result
        filtered.to_csv(output.filtered_list, header=True, index=False, sep='\t')
        out_data = filtered.loc[:,["name", "downsample"]]
        normalizing_sample = out_data[out_data.downsample == 1] # which one is hte normalizing sample
        print(normalizing_sample)
        out_data = out_data[out_data.downsample < 1] #only keep samples that need to be downsampled and skip the ones that wont be
        out_data["downsample"] = (round(out_data["downsample"], 8)*100000000).astype(int) #round the sample to the nearest number, I will add the . later
        out_data.to_csv(output.downsample_list, header=False, index=False, sep='\t') #write everything out to a csv

#ATAC sequencing involves a transposase that slightly shifts the reads when 
#it is used. In order to correct this we need to 5, 4 offset the reads so that they line up with the actual accessibility locations
rule shiftReads:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
        index = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
    conda:
        "../envs/deeptools.yaml"
    threads: 4 
    message:
        """--- shifting reads---"""
    shell:
        """
        alignmentSieve --bam {input.bamfile} --outFile {output.bamfile} --numberOfProcessors {threads} --ATACshift

        """

#many commands require that there is index for the bamfile that you are using. TODO: this command should
#eventually be pulled together into the command that needs this index which is hte shifting rule  above
rule index:
    input:
        deduplicated = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
        full_bam = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
    output:
        indexed = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
        indexed_full = sample_work_path + "/bamfiles/{merged_sample}_merged.bam.bai",
    conda:
        "../envs/sambamba.yaml"
    threads: 4 
    message:
        """--- indexing reads---"""
    shell:
        """
        sambamba index -t {threads} {input.full_bam} {output.indexed_full}
        sambamba index -t {threads} {input.deduplicated} {output.indexed}
        
        """

#Filter out reads that match the quality flags that we want in samtools
rule samtoolsQuality:
     input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam",
     output:
        #TODO add in temp()
        quality = temp(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam"),
     params:
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv"
     conda:
        "../envs/samtools_env.yaml"
     threads: 4
     message:
        """--- deduplicating reads---"""
     shell:
        """
        samtools view -h -@ 4 -b -F 1804 {input.bamfile} > {output.quality}
        samtools view -@ 4 -c {output.quality}| xargs -I{{}} echo "{wildcards.merged_sample},satisfy_quality,"{{}}$'\n'  >> {params.readsfile}
        """
        # Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
        # Retain properly paired reads -f 2

#deduplicate the reads that we have. This will get rid of potential PCR duplicates
#After this step you might have a large difference in reads if some samples have a high percentage of duplicates 
rule dedup:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
    output:
        #TODO add in temp()
        deduplicated = temp(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam"),
    params:
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv",
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
        """--- deduplicating reads---"""
    shell:
         """
         sambamba markdup -t {threads} {input.bamfile} {output.deduplicated}
         samtools view -@ 4 -F 1024 -c {output.deduplicated}| xargs -I{{}} echo "{wildcards.merged_sample},no_duplicates,"{{}}$'\n'  >> {params.readsfile}
         """

#Remove the mitochondrial reads This requires that the mitochondrial reads have the  name chrM so that you can grep them with an inverted search
#This works for our current aligner and references, but if you get zero mitochondrial reads filtered out its possible that the mitochondrial chromosome
#has a name that is different from "chrM"
rule removeMitochondrial:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_sorted.bam", #change this to _merged if merging lanes
    output:
        #TODO add in temp()
        bamfile = temp(sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam"),
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv"
    conda:
        "../envs/samtools_env.yaml"
    message:
        """--- Removing mitochondrial reads---"""
    log:
        totalreads = sample_work_path + "/bamfiles/{merged_sample}_readswmito.txt"
    threads: 4 
    shell:
        """
        samtools view -@ 4 -c {input.bamfile}| xargs -I{{}} echo "{wildcards.merged_sample},starting_reads,"{{}}$'\n' >> {output.readsfile}
        samtools view -h -@ {threads} {input.bamfile}| grep -v chrM| samtools sort -@ {threads} -O BAM > {output.bamfile}
        samtools view -@ 4 -c {output.bamfile}| xargs -I{{}} echo "{wildcards.merged_sample},no_ChrM,"{{}}$'\n' >> {output.readsfile}

        """

#WE should be able to remove this part for anything that is not run on mutliple lanes, however not entirely sure
#Will have to check once the data comes in
#This allows us to merge samples from multiple lanes if necessary and takes the number of lanes from the
#configuration file.
rule mergeBam:
    input:
        lane1 = sample_work_path + "/bamfiles/{merged_sample}" + LANES[0] + "_sorted.bam",
        lane2 = sample_work_path + "/bamfiles/{merged_sample}" + LANES[1] + "_sorted.bam",
    output:
        merged_bamfile = temp(sample_work_path + "/bamfiles/{merged_sample}_merged.bam"),
    conda:
        "../envs/samtools_env.yaml"
    threads: 8 
    message:
         """--- sorting the bamfile ---"""
    shell:
         """ cp {input.lane1} {output.merged_bamfile} """
# should not be double counting, however if lanes are different can use this command
##         """samtools merge -@ {threads} {output.merged_bamfile} {input.lane1} {input.lane2}"""



