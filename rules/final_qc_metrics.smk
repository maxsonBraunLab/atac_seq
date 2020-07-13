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


