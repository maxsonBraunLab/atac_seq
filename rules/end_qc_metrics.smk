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


rule peak_catalog_counts_no_downsmpl:
    input:
        read_catalog = sample_work_path + "/bamfiles/reads_catalog_nodownsample.bed",
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
    output:
        read_count = sample_work_path + "/fully_filtered/{merged_sample}_read_catalog_nodownsample_counts.bed",
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam",
    params:
        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + "{merged_sample}",
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

rule peak_catalog_no_downsmpl:
    input:
        peaks = expand(sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES),
    output:
        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_nodownsample.bed",
    params:
        merged_bed = sample_work_path + "/bamfiles/all_broad_peaks_nodownsample.bed",
        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_intervals_nodownsample.bed",
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






#    params:
#        merged_bed = sample_work_path + "/bamfiles/all_broad_peaks_nodownsample.bed",
#        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_nodownsample.bed",
#        present_in_number = 2,
#        blacklist = blacklist_file,
#        genome = genome_name,
#        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + '\t'.join(sorted(expand("{sample_id}", sample_id=MERGED_SAMPLES))),
#        peaks_input = " ".join(sorted(expand(sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES))),
#        script_file = "./scripts/computing_peak_catalog.R",
#        bams_input = " ".join(sorted(expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam", merged_sample=MERGED_SAMPLES))),
#        #we are sorting both of our bamfiles and our header. This will order the headers as well as the bams in a way that makes them group together so that it is easier for analysis
#
#    conda:
#        "../envs/peaks_catalog.yaml"
#    shell:
#        """
#        echo "" > {params.merged_bed}
#        cat {params.peaks_input} >> {params.merged_bed}
#        Rscript --vanilla {params.script_file} {params.merged_bed} {params.blacklist} {params.present_in_number} {params.reads_catalog_bed} {params.genome}
#        echo {params.header}
#        echo {params.bams_input}
#        """

#bedtools multicov
#sed "1i {params.header}" {output.reads_catalog} > {output.reads_catalog}



rule MACS_no_downsmpl:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam"
    output:
        summits = sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample/{merged_sample}_macs_peaks.broadPeak"
    conda:
        "../envs/MACS.yaml"
    threads: 4
    params:
        genome = macs_genome_size, 
        name = "{merged_sample}_macs",
        directory = sample_work_path + "/bamfiles/{merged_sample}_macsout_nodownsample"
    message:
        """--- running MACS ---"""
    shell:
        """
        macs2 callpeak --treatment {input.bamfile} --name {params.name} --format BAMPE --gsize {params.genome} --outdir {params.directory} --broad
        """


