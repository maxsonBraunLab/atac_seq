rule peak_catalog:
    input:
        peaks = expand(sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES),
        bamfiles = expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam", merged_sample=MERGED_SAMPLES),
    output:
        reads_catalog = sample_work_path + "/bamfiles/reads_catalog_counts.bed",
    params:
        merged_bed = sample_work_path + "/bamfiles/all_broad_peaks.bed",
        reads_catalog_bed = sample_work_path + "/bamfiles/reads_catalog_intervals.bed",
        present_in_number = 2,
        blacklist = blacklist_file,
        header = ", , , ," + ", ".join(expand("{sample_id}", sample_id=MERGED_SAMPLES)) + '\n',
        peaks_input = " ".join(expand(sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak", merged_sample=MERGED_SAMPLES)),
        script_file = "./scripts/computing_peak_catalog.R",
        bams_input = " ".join(expand(sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam", merged_sample=MERGED_SAMPLES)),
    conda:
        "../envs/peaks_catalog.yaml"
    shell:
        """
        echo "" > {params.merged_bed}
        cat {params.peaks_input} >> {params.merged_bed}
        Rscript --vanilla {params.script_file} {params.merged_bed} {params.blacklist} {params.present_in_number} {params.reads_catalog_bed}
        bedtools multicov -bams {params.bams_input} -bed {params.reads_catalog_bed} >> {output.reads_catalog}
        """

#bedtools multicov
#sed "1i {params.header}" {output.reads_catalog} > {output.reads_catalog}


#rule HOMER:
#    input:
#        summits = sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.narrowPeak"
#    output:
#        all_motifs = sample_work_path + "/bamfiles/{merged_sample}_homerout/homerMotifs.all.motifs"
#    conda:
#        "../envs/HOMER.yaml"
#    threads: 4
#    params:
#        genome = "/home/groups/MaxsonLab/indices/mm10/mm10.fa",
#        size = "200",
#        directory = sample_work_path + "/bamfiles/{merged_sample}_homerout/"
#    message:
#        """--- running HOMER---"""
#    shell:
#        """
#        findMotifsGenome.pl {input.summits} {params.genome} {params.directory} -size {params.size} -p {threads}
#        """



rule MACS:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam"
    output:
        summits = sample_work_path + "/bamfiles/{merged_sample}_macsout/{merged_sample}_macs_peaks.broadPeak"
    conda:
        "../envs/MACS.yaml"
    threads: 4
    params:
        genome = "1.87e9",
        name = "{merged_sample}_macs",
        directory = sample_work_path + "/bamfiles/{merged_sample}_macsout"
    message:
        """--- running MACS ---"""
    shell:
        """
        macs2 callpeak --treatment {input.bamfile} --name {params.name} --format BAMPE --gsize {params.genome} --outdir {params.directory} --broad
        """


