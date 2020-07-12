rule non_downsampled_bedGraph:
    input:
        non_downsampled = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam",
    output:
        non_downsampled_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_non_downsampled.bedGraph",
    conda:
        "../envs/peaks_catalog.yaml"
    message:
        """--- making non_downsampled bedGraph---"""
    shell:
        """
        bedtools genomecov -ibam {input.non_downsampled} -bg > {output.non_downsampled_bedGraph}

        """

rule non_downsampled_sort:
    input:
        non_downsampled_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_non_downsampled.bedGraph",
    output:
        non_downsampled_sorted_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_non_downsampled_sorted.bedGraph",
    conda:
        "../envs/bedSort.yaml"
    message:
        """--- sorting non_downsampled bedGraph---"""
    shell:
        """
        bedSort {input.non_downsampled_bedGraph} {output.non_downsampled_sorted_bedGraph}

        """

rule non_downsampled_bedGraphToBigWig:
    input:
        non_downsampled_sorted_bedGraph = sample_work_path + "/bamfiles/{merged_sample}_non_downsampled_sorted.bedGraph",
    output:
        non_downsampled_bigwig = sample_work_path + "/non_downsampled_bigwigs/{merged_sample}_non_downsampled.bw",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    params:
        chrom_size = "/home/groups/MaxsonLab/indices/mm10/mm10.chrom.sizes"
    message:
        """--- making non_downsampled bigwig ---"""
    shell:
        """
        bedGraphToBigWig {input.non_downsampled_sorted_bedGraph} {params.chrom_size} {output.non_downsampled_bigwig}
       
        """

