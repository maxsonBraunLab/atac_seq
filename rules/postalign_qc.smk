rule makeBigwig:
    input:
        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
        deduplicated = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
    output:
        bigwig = sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_smooth.bw",
        bigwig_rough = sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_rough.bw",

    conda:
        "../envs/deeptools.yaml"
    threads: 8 
    message:
        """--- making bigwig---"""
    shell:
        """
        bamCoverage --bam {input.quality} -o {output.bigwig} --numberOfProcessors {threads} --skipNAs --binSize 5 --smoothLength 15  --normalizeUsing RPKM
        bamCoverage --bam {input.quality} -o {output.bigwig_rough} --numberOfProcessors {threads} --skipNAs --binSize 5 --normalizeUsing RPKM

        """
#rule bamcoverage bedtools multicov [OPTIONS] -bams BAM1 BAM2 BAM3 ... BAMn -bed  <BED/GFF/VCF>
# This uses the r_analysis environment
#

rule downSample:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
        downsample_list = sample_work_path + "/bamfiles/downsample_list.txt"
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReadsDownSample.bam",
    conda:
        "../envs/sambamba.yaml"
    threads: 4 
    message:
        """---normalizing reads---"""
    shell:
        """
        export PATH=$PATH:/opt/installed/samtools-1.6/bin/
        samtools view -b -h -s 1.$(grep {wildcards.merged_sample} {params.subsample_file}| awk '{{print $2}}') {input.bamfile} {output.bamfile}        
        
        """

rule createDownsample:
    input:
        bamfile = expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam", merged_sample=MERGED_SAMPLES),
        readsfile = expand(sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv", merged_sample=MERGED_SAMPLES),
        index = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam.bai",
    output:
        downsample_list = sample_work_path + "/bamfiles/downsample_list.txt"
    message:
        """--calculating normalization factor to subsample by---"""
    run:
        import pandas as pd
        dataframes = []
        for file in input.readsfile:
            dataframes.append(pd.read_csv(file, names=["name","description","reads"]))
        merged_dataframe = pd.concat(dataframes, ignore_index=True)
        filtered = merged_dataframe[merged_dataframe.description=="satisfy_quality"]
        filtered['downsample'] = filtered['reads'].min() / filtered['reads']
        out_data = filtered.loc[:,["name", "downsample"]]
        out_data.to_csv(output.downsample_list, header=False, index=False, sep='\t')


rule shiftReads:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
        indexed = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
        index = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam.bai",
    conda:
        "../envs/deeptools.yaml"
    threads: 4 
    message:
        """--- shifting reads---"""
    shell:
        """
        alignmentSieve --bam {input.bamfile} --outFile {output.bamfile} --numberOfProcessors {threads} --ATACshift

        export PATH=$PATH:/opt/installed/samtools-1.6/bin/
        samtools index -b {output.bamfile} {output.index}        

        """

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

#rule countReads:
#     input:
#        merged = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
#        dedup = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
#        no_mito = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam",
#        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam"
#     output:
#        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv"
#     conda:
#        "../envs/samtools_env.yaml"
#     threads: 4
#     message:
#        """--- counting reads ---"""
#     shell:
#        """
#        samtools view -@ 4 -c {input.merged}| xargs -I{{}} echo "{wildcards.merged_sample},starting_reads,"{{}}$'\n' >> {output.readsfile}
#        samtools view -@ 4 -c {input.no_mito}| xargs -I{{}} echo "{wildcards.merged_sample},no_ChrM,"{{}}$'\n' >> {output.readsfile}
#        samtools view -@ 4 -c {input.dedup}| xargs -I{{}} echo "{wildcards.merged_sample},no_duplicates,"{{}}$'\n'  >> {output.readsfile}
#        samtools view -@ 4 -c {input.quality}| xargs -I{{}} echo "{wildcards.merged_sample},satisfy_quality,"{{}}$'\n'  >> {output.readsfile}
#        """


rule samtoolsQuality:
     input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam",
     output:
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
 
rule dedup:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
    output:
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

rule removeMitochondrial:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_sorted.bam", #change this to _merged if merging lanes
    output:
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
         """samtools merge -@ {threads} {output.merged_bamfile} {input.lane1} {input.lane2}"""



