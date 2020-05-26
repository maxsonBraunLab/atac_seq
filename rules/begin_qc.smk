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

rule sortbam:
    input:
        bamfile = sample_work_path + "/bamfiles/{sample}.bam",
    output:
        sorted_bamfile = sample_work_path + "/bamfiles/{sample}_sorted.bam",
    conda:
        "../envs/samtools_env.yaml"
    params: outdirectory = sample_work_path + "/bamfiles" #this is where temporary compression files are put
    threads: 8 
    message:
        """--- sorting the bamfile ---"""
    shell:
        """samtools sort -T {params.outdirectory} -@ {threads} -o {output.sorted_bamfile} {input.bamfile} """


rule sam2bam:
    input:
        samfile = sample_work_path + "/samfiles/{sample}.sam",
    output:
        bamfile = sample_work_path + "/bamfiles/{sample}.bam",
    conda:
        "../envs/samtools_env.yaml"
    threads: 4 
    message:
        """--- converting {input.samfile} to {output.bamfile} ---"""
    shell:
        """samtools view -bh -@ {threads} -o {output.bamfile} {input.samfile}"""

rule align:
    input:
        forward_paired = sample_work_path + "/trimmed/{sample}_1_fastqc_tpaired.fastq",
        reverse_paired = sample_work_path + "/trimmed/{sample}_2_fastqc_tpaired.fastq",

    output:
        samfile = temp(sample_work_path + "/samfiles/{sample}.sam"),
    params:
        reference = reference_align_genome
    conda:
        "../envs/bwa.yaml"
    message:
        """--- Aligning reads to reference."""
    threads: 8 
    shell:
        """bwa mem -t {threads} {params.reference} {input.forward_paired} {input.reverse_paired} > {output.samfile} """


rule trimming:
    input:
        forward = raw_sample_folder + prefix + "{sample}" + readforward_postfix,
        reverse = raw_sample_folder + prefix + "{sample}" + readreverse_postfix,
    output:
        forward_paired = temp(sample_work_path + "/trimmed/{sample}_1_fastqc_tpaired.fastq"),
        reverse_paired = temp(sample_work_path + "/trimmed/{sample}_2_fastqc_tpaired.fastq"),
        forward_unpaired = temp(sample_work_path + "/trimmed/{sample}_1_fastqc_tunpaired.fastq"),
        reverse_unpaired = temp(sample_work_path + "/trimmed/{sample}_2_fastqc_tunpaired.fastq"),
    params:
        adapter = adapters_file #this variable is defined in the Snakefile
    conda:
        "../envs/trim.yaml"
    threads: 4
    message:
        """--- Trimming."""
    shell:
        """trimmomatic PE {input.forward} {input.reverse} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired}   ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"""


rule fastqc:
    input:
        forward = raw_sample_folder + prefix + "{sample}" + readforward_postfix,
        reverse = raw_sample_folder + prefix + "{sample}" + readreverse_postfix,
    output:
        forward = sample_work_path + "/fastqc/{sample}/{sample}_1_fastqc.zip",
        reverse = sample_work_path + "/fastqc/{sample}/{sample}_2_fastqc.zip",
    params: outdir= sample_work_path + "/fastqc/{sample}",
    conda:
        "../envs/fastqc.yaml"
    threads: 4
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """
        fastqc --outdir  {params.outdir} --extract  -f fastq {input.forward} {input.reverse}
        """


