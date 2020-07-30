rule merge_catalog_counts:
    input:
         read_count = expand(sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_counts.bed", merged_sample=MERGED_SAMPLES),
    output:
        read_count = sample_work_path + "/fully_filtered/all_read_catalog_counts.bed",
    run:
        import pandas as pd
        dataframes = []
        for file in input.read_count:
            dataframes.append(pd.read_csv(file, sep='\t'))
        merged = pd.concat(dataframes, axis=1)
        merged = merged.loc[:,~merged.columns.duplicated()]   
        merged = merged.sort_index(axis=1)
        merged.to_csv(output.read_count, header=True, index=False, sep='\t')
 

rule peak_catalog_counts:
    input:
        read_catalog = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam",
    output:
        read_count = sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_counts.bed",
    params:
        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + "{merged_sample}",
    conda:
        "../envs/peaks_catalog.yaml"
    shell:
        """
        bedtools multicov -bams {input.bamfile} -bed {input.read_catalog} >> {output.read_count}
        sed -i "1i {params.header}" {output.read_count}
       
        """

rule annotate_peak_names:
    input:
        bedfile = sample_work_path + "/fully_filtered/reads_catalog_intervals.bed",
    output:
        read_count = sample_work_path + "/fully_filtered/{merged_sample}_peaks_information.txt",
    params:
        genome = genome_name,
    conda:
        "../envs/HOMER.yaml"
    shell:
        """
        annotatePeaks.pl {input.bedfile}  {params.genome} > {output.read_count}
        """


rule peak_catalog:
    input:
        peaks = expand(sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES),
        #bamfiles = expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam", merged_sample=MERGED_SAMPLES),
    output:
        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
        sample_condition_catalog = expand(sample_work_path + "/bamfiles/reads_catalog_intervals.bed_{sample_condition}.bed",sample_condition=SAMPLE_CONDITIONS),

    params:
        merged_bed = sample_work_path + "/bamfiles/all_broad_peaks.bed",
        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
        present_in_number = 2,
        blacklist = blacklist_file,
        genome = genome_name,
        metadata = metadata_file,
        regex = remove_filepath_regex,

        peaks_input = " ".join(sorted(expand(sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES))),
        script_file = "./scripts/computing_peak_catalog_by_replicate.R",

    conda:
        "../envs/peaks_catalog.yaml"
    shell:
        """
        echo "" > {params.merged_bed}
        cat {params.peaks_input} >> {params.merged_bed}
        Rscript --vanilla {params.script_file} {params.merged_bed} {params.blacklist} {params.present_in_number} {params.reads_catalog_bed} {params.genome} {params.metadata} {params.regex} 
        """


rule MACS:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_downSample.bam"
    output:
        summits = sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak"
    conda:
        "../envs/MACS.yaml"
    threads: 4
    params:
        genome = macs_genome_size, 
        name = "{merged_sample}_macs",
        directory = sample_work_path + "/bamfiles/{merged_sample}_macsout"
    message:
        """--- running MACS ---"""
    shell:
        """
        macs2 callpeak --treatment {input.bamfile} --name {params.name} --format BAMPE --gsize {params.genome} --outdir {params.directory} --broad
        """


