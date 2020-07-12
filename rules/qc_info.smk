rule atac_qc:
    input:
        merged_bamfile = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
    output:

    conda:
        "../envs/macs.yaml"
    threads: 4 
    script: "../scripts/atac_qc.R"
    message:
        """---make qc plots for individual samples ---"""
    shell:
        """  """


