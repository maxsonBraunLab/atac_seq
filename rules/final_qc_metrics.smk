#move all the end results files to a results folder
rule move_files:
    input:
        sample_informtion_table = sample_work_path + "/results/sample_quality_information.tsv"
        correlation_plot = sample_work_path + "/fully_filtered/correllation_plot.png",
        intersected_read_catalog = sample_work_path + "/bamfiles/intersected_catalog.bed",
        read_catalog = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
        read_count = sample_work_path + "/fully_filtered/all_read_catalog_nodownsample_counts_ondownsamplepeaks.bed",
    params:
        result_dir = sample_work_path + "/results"
    output:
        sample_informtion_table = sample_work_path + "/results/sample_quality_information.tsv"
        correlation_plot = sample_work_path + "/results/correlation_plot.png",
        intersected_read_catalog = sample_work_path + "/results/intersected_catalog.bed",
        read_catalog = sample_work_path + "/results/reads_catalog_intervals.bed",
        read_count = sample_work_path + "/results/all_read_catalog_nodownsample_counts_ondownsamplepeaks.bed",
    shell:
        """
        mkdir {params.result_dir}
        cp {input.sample_information_table} {output.sample_information_table}
        cp {input.correlation_plot} {output.correlation_plot}
        cp {input.intersected_read_catalog} {output.intersected_read_catalog}
        cp {input.read_catalog} {output.read_catalog}
        cp {input.read_count} {output.read_count}
        """

#create a table that has quality information for each one of the samples
rule create_quality_table:
    input:
        read_depths_qc = expand(sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv", merged_sample=MERGED_SAMPLES),
        done = sample_work_path + "/fully_filtered/{merged_sample}_count_peaks_done.txt",
    params:
        fully_filtered_path = sample_work_path + "/fully_filtered/*.csv",
        sample_information_all = sample_work_path + "/fully_filtered/all_sample_information_nospaces.csv",
    output:
        sample_informtion_table = sample_work_path + "/results/sample_quality_information.tsv",
    shell:
        """
        cat {params.fully_filtered_path} | sed '/^\s*$/d' > {params.sample_information_all}
        Rscript --vanilla {params.sample_information_all} {output.sample_information_table}  
        """


rule plotCorrelation:
    input:
        bigwig_rough = expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_rough.bw", merged_sample=MERGED_SAMPLES),
    params:
        bigwigs = " ".join(sorted(expand(sample_work_path + "/fully_filtered/{merged_sample}_tracks_5window_rough.bw", merged_sample=MERGED_SAMPLES))),
        names = " ".join(sorted(expand("{merged_sample}", merged_sample=MERGED_SAMPLES))),
    output:
        scores = sample_work_path + "/fully_filtered/sample_scores.npz", 
        correlation_plot = sample_work_path + "/fully_filtered/correllation_plot.png",
    conda:
        "../envs/deeptools.yaml"
    threads: 8 
    message:
        """Creating multibigwig"""
    shell:
        """
        multiBigwigSummary bins -b {params.bigwigs} --labels {params.names} -out {output.scores} -p {threads}
        
        plotCorrelation -in {output.scores} --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts"  --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.correlation_plot} 
        """

#This runs the bedtools multi intersect command this command tells us which samples contribute to each of the
#intervals, this allows us to tell which intervals are only present in certain sample conditions
rule multi_intersect:
    input:
        expand(sample_work_path + "/bamfiles/reads_catalog_intervals.bed_{sample_condition}_formatted.bed", sample_condition=SAMPLE_CONDITIONS),
    output:
        bedfile = sample_work_path + "/bamfiles/intersected_catalog.bed",
    params:
        bedfiles =  " ".join( sorted(expand(sample_work_path + "/bamfiles/reads_catalog_intervals.bed_{sample_condition}_formatted.bed", sample_condition=sorted(SAMPLE_CONDITIONS)))),
        names = " ".join( sorted(SAMPLE_CONDITIONS) ),
    shell:
        """
        export PATH=$PATH:/home/groups/MaxsonLab/software/bedtools/
        bedtools multiinter -i {params.bedfiles} -header -names {params.names} > {output.bedfile}

        """

#This runs a trick with bedtools intersect that allows us to determine
#which of the overall intervals an interval in a sample is from and returns that interval
#If we don't do this bedtools multiintersect will not provide contiguous intervals
rule get_union_interval_from_sample_interval:
    input:
        sample_condition_catalog = sample_work_path + "/bamfiles/reads_catalog_intervals.bed_{sample_condition}.bed",
        read_catalog = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
    output:
        sample_condition_catalog = sample_work_path + "/bamfiles/reads_catalog_intervals.bed_{sample_condition}_formatted.bed",
    shell:
        """
        export PATH=$PATH:/home/groups/MaxsonLab/software/bedtools/
        bedtools intersect -wa -a {input.read_catalog} -b {input.sample_condition_catalog} > {output.sample_condition_catalog}

        """

#This command counts the number of MACS peaks from the non downsampled bam files for usage in the 
#quality catalog that we will be making
rule count_macs_peaks:
    input:
        done = sample_work_path + "/fully_filtered/{merged_sample}_information_done.txt",
        summits = sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample/{merged_sample}_macs_peaks.broadPeak",
    output:
        done = sample_work_path + "/fully_filtered/{merged_sample}_count_peaks_done.txt",
    params:
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv",
    shell:
        """
        wc -l {input.summits}| awk '{{print $1}}'| xargs -I {{}} echo "{wildcards.merged_sample},unfiltered_peaks,"{{}}$'\n' >> {params.readsfile} 
        touch {output.done}
        """

#This rule provides a large amount of information that can be used for the quality information catalog
#This information is all stored in a sample specific file in the path sample_work_path/fully_filtered/{merged_sample}_read_depths.csv
#You can create a concatenation of all of these files by going into the fully filtered folder and running cat */*.csv 
rule get_info:
    input:
        peaks = sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample/{merged_sample}_macs_peaks.broadPeak",
        bamfile = sample_work_path + "/bamfiles/{merged_sample}.bam",
        read_count = sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_nodownsample_counts_ondownsamplepeaks.bed",
        #read_count = sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_nodownsample_counts.bed",
        filtered_bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
    output:
        done = sample_work_path + "/fully_filtered/{merged_sample}_information_done.txt",
        flagstat = sample_work_path + "/fully_filtered/{merged_sample}_flagstat.txt",
    params:
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv",
        flagstat = sample_work_path + "/fully_filtered/{merged_sample}_flagstat.txt",

    shell:
        """
        export PATH=$PATH:/opt/installed/samtools-1.6/bin/
        samtools flagstat {input.bamfile} > {params.flagstat}
        grep total {params.flagstat}|awk '{{print $1}}'|xargs -I {{}} echo "{wildcards.merged_sample},total_reads_trimd,"{{}}$'\n' >> {params.readsfile} 
        grep "[0-9] mapped" {params.flagstat} |awk '{{print $1}}' | xargs -I {{}} echo "{wildcards.merged_sample},mapping_reads,"{{}}$'\n' >> {params.readsfile} 
        awk '{{sum +=$7}} END{{print sum}}' {input.read_count}| xargs -I {{}} echo "{wildcards.merged_sample},reads_in_peaks,"{{}}$'\n' >> {params.readsfile} 
        touch {output.done}
        """  


