

#rule annotate_peak_names_no_downsmpl:
#    input:
#        bedfile = sample_work_path + "/bamfiles/intersected_bedfile_nodownsample.bed",
#    output:
#        peak_information_catalog = sample_work_path + "/fully_filtered/intersected_bedfile_nodownsample_peakinfo.txt",
#    params:
#        genome = genome_name,
#    conda:
#        "../envs/HOMER.yaml"
#    shell:
#        """
#        annotatePeaks.pl {input.bedfile} {params.genome} > {output.peak_information_catalog}
#        """

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
        read_count = sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_nodownsample_counts.bed",
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


