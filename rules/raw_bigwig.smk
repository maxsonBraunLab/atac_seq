rule raw_bedGraph:
    input:
        raw_bam = sample_work_path + "/bamfiles/{merged_sample}_sorted.bam",
    output:
        raw_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_raw.bedGraph",
    conda:
        "../envs/peaks_catalog.yaml"
    message:
        """--- making raw bedGraph---"""
    shell:
        """
        bedtools genomecov -ibam {input.raw_bam} -bg > {output.raw_bedGraph}

        """

rule raw_sort:
    input:
        raw_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_raw.bedGraph",
    output:
        raw_sorted_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_raw_sorted.bedGraph",
    conda:
        "../envs/bedSort.yaml"
    message:
        """--- sorting raw bedGraph---"""
    shell:
        """
        bedSort {input.raw_bedGraph} {output.raw_sorted_bedGraph}

        """

rule raw_bedGraphToBigWig:
    input:
        raw_sorted_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_raw_sorted.bedGraph",
    output:
        raw_bigwig = sample_work_path + "/raw_bigwigs/{merged_sample}_raw.bw",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    params:
        chrom_size = "/home/groups/MaxsonLab/indices/mm10/mm10.chrom.sizes"
    message:
        """--- making raw bigwig ---"""
    shell:
        """
        bedGraphToBigWig {input.raw_sorted_bedGraph} {params.chrom_size} {output.raw_bigwig}
       
        """

