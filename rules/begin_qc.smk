#rule fastqscreen:
#    input:
#        "../samples/trimmed/{sample}_t.fastq"
#    output:
#        "../samples/fastqscreen/{sample}/{sample}_t_screen.html",
#        "../samples/fastqscreen/{sample}/{sample}_t_screen.png",
#        "../samples/fastqscreen/{sample}/{sample}_t_screen.txt"
#    params:
#        conf = config["conf"]
#    conda:
#        "../envs/fastqscreen.yaml"
#    shell:
#        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input}"""
#
#rule shiftReads:
#    input:
#        bamfile = "../samples/bamfiles/{sample}_rmChrM_dedup.bam",
#    output:
#        deduplicated = "../samples/bamfiles/{sample}_rmChrM_dedup_shiftedReads.bam"
#    params:
#        #reference = config["reference_align_genome"], 
#    conda:
#        "../envs/atacseqqc.yaml"
#    message:
#        """--- shifting reads---"""
#    shell:
#        """ """
#
#rule dedup:
#    input:
#        bamfile = "../samples/bamfiles/{sample}_rmChrM.bam",
#    output:
#        deduplicated = "../samples/bamfiles/{sample}_rmChrM_dedup.bam"
#    params:
#        #reference = config["reference_align_genome"], 
#    conda:
#        "../envs/picard.yaml"
#    message:
#        """--- deduplicating reads---"""
#    shell:
#        """ """
#
#rule removeMitochondrial:
#    input:
#        bamfile = "../samples/bamfiles/{sample}.bam",
#    output:
#        bamfile = "../samples/bamfiles/{sample}_rmChrM.bam",
#    params:
#        #reference = config["reference_align_genome"], 
#    conda:
#        "../envs/samtools.yaml"
#    message:
#        """--- Removing mitochondrial reads---"""
#    shell:
#        """ """
#
#rule sortbam:
#    input:
#        bamfile = "../samples/bamfiles/{sample}.bam",
#    output:
#        bamfile = "../samples/bamfiles/{sample}.bam",
#    params:
#        #reference = config["reference_align_genome"], 
#    conda:
#        "../envs/samtools.yaml"
#    message:
#        """--- indexing the bamfile ---"""
#    shell:
#        """ """
#
#rule sam2bam:
#    input:
#        samfile = "../samples/samfiles/{sample}.sam",
#    output:
#        bamfile = "../samples/bamfiles/{sample}.bam",
#    params:
#        #reference = config["reference_align_genome"], 
#    conda:
#        "../envs/samtools.yaml"
#    message:
#        """--- converting {input.samfile} to {output.bamfile} ---"""
#    shell:
#        """ """

rule align:
    input:
        forward_paired = "../samples/trimmed/{sample}_1_fastqc_tpaired.fastq",
        reverse_paired = "../samples/trimmed/{sample}_2_fastqc_tpaired.fastq",

    output:
        samfile = "../samples/samfiles/{sample}.sam",
    params:
        #reference = config["reference_align_genome"], 
        reference = "/home/groups/MaxsonLab/indices/GRch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
    conda:
        "../envs/bwa.yaml"
    message:
        """--- Aligning reads to reference."""
    threads: 8 
    #takes a look at how many cores are available to snakemake and give 75 % of them to this rule as it is the slowest one
    shell:
        """bwa mem -t 8 {params.reference} {input.forward_paired} {input.reverse_paired} > {output.samfile} """


rule trimming:
    input:
        forward = "../samples/raw/{sample}_1.fastq",
        reverse = "../samples/raw/{sample}_2.fastq",

    output:
        forward_paired = "../samples/trimmed/{sample}_1_fastqc_tpaired.fastq",
        reverse_paired = "../samples/trimmed/{sample}_2_fastqc_tpaired.fastq",
        forward_unpaired = "../samples/trimmed/{sample}_1_fastqc_tunpaired.fastq",
        reverse_unpaired = "../samples/trimmed/{sample}_2_fastqc_tunpaired.fastq",
    params:
        adapter = config["adapters_file"]
    conda:
        "../envs/trim.yaml"
    threads: 4
    message:
        """--- Trimming."""
    shell:
        """trimmomatic PE {input.forward} {input.reverse} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired}   ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"""


rule fastqc:
    input:
        forward = "../samples/raw/{sample}_1.fastq",
        reverse = "../samples/raw/{sample}_2.fastq",
    output:
        forward = "../samples/fastqc/{sample}/{sample}_1_fastqc.zip",
        reverse = "../samples/fastqc/{sample}/{sample}_2_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    threads: 4
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """
        fastqc --outdir  ../samples/fastqc/{wildcards.sample} --extract  -f fastq {input.forward} {input.reverse}
        """


