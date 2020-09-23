#This rule merges all of the catalog counts from the individual samples into one catalog for all of the samples
rule merge_catalog_counts_no_downsmpl:
    input:
        read_count = expand("samples/fully_filtered/{sample}_read_catalog_nodownsample_counts.bed", sample=SAMPLE),
    output:
        read_count = "samples/fully_filtered/all_read_catalog_nodownsample_counts.bed",
    run:
        import pandas as pd
        dataframes = []
        for file in input.read_count:
            dataframes.append(pd.read_csv(file, sep='\t'))
        merged = pd.concat(dataframes, axis=1)
        merged = merged.loc[:,~merged.columns.duplicated()]   
        merged = merged.sort_index(axis=1)
        merged.to_csv(output.read_count, header=True, index=False, sep='\t')


#This creates the peak catalog for an individual sample by counting the coverage individually
#When counting reads in intervals all together bamcov had issues and would drop low count samples
#This allows us to include low count samples
rule peak_catalog_counts_no_downsmpl:
    input:
        read_catalog = "samples/bamfiles/reads_catalog_intervals_nodownsample.bed",
        bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads.bam",
    output:
        read_count = "samples/fully_filtered/{sample}_read_catalog_nodownsample_counts.bed",
        bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam",
    params:
        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + "{sample}",
    conda:
        "../envs/peaks_catalog.yaml"
    threads: 8 
    shell:
        """
        export PATH=$PATH:/opt/installed/samtools-1.6/bin/
        samtools sort -@ 4 {input.bamfile} > {output.bamfile}
        samtools index {output.bamfile}
         
        bedtools multicov -bams {output.bamfile} -bed {input.read_catalog} >> {output.read_count}
        sed -i "1i {params.header}" {output.read_count}
       
        """

#This rule has two outputs, one is the overall peak catalog and another is the peak catalog for each sample condition
#This can later be input into the bedtools commands to allow us to assign peaks in the sample condition catalogs to the 
#corresponding peak in the overall peak catalog so that we can run bedtools multiintersect
rule peak_catalog_no_downsmpl:
    input:
        peaks = expand("samples/bamfiles/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak", sample=SAMPLE),
    output:
        reads_catalog_bed = "samples/bamfiles/reads_catalog_intervals_nodownsample.bed",
    params:
        merged_bed = "samples/bamfiles/all_broad_peaks_nodownsample.bed",
        present_in_number = 2,
        blacklist = config["blacklist_file"],
        genome = config["genome_name"],
        metadata = config["metadata_file"],
        regex = "", #This assumes that the sample names will exactly match 
        #What is on the metadata file. This feature should be removed as
        #its implementation is lacking

        peaks_input = " ".join(sorted(expand("samples/bamfiles/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak", sample=SAMPLE))),
        script_file = "./scripts/computing_peak_catalog_by_replicate.R",

    conda:
        "../envs/peaks_catalog.yaml"
    shell:
        """
        echo "" > {params.merged_bed}
        cat {params.peaks_input} >> {params.merged_bed}
        Rscript --vanilla {params.script_file} {params.merged_bed} {params.blacklist} {params.present_in_number} {output.reads_catalog_bed} {params.genome} {params.metadata} {params.regex} 
        """

#This command runs MACS on the non downsampled bam files
rule MACS_no_downsmpl:
    input:
        bamfile = "samples/bamfiles/{sample}_rmChrM_dedup_quality_shiftedReads.bam"
    output:
        summits = "samples/bamfiles/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak"
    conda:
        "../envs/MACS.yaml"
    threads: 4
    params:
        genome = config["macs_genome_size"], 
        name = "{sample}_macs",
        directory = "samples/bamfiles/{sample}_macsout_nodownsample"
    message:
        """--- running MACS ---"""
    shell:
        """
        macs2 callpeak --treatment {input.bamfile} --name {params.name} --format BAMPE --gsize {params.genome} --outdir {params.directory} --broad
        """


