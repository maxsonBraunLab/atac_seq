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
#error: reverse is reserved for internal use
    input:
        forward_strand = raw_sample_folder + prefix + "{sample}" + readforward_postfix,
        reverse_strand = raw_sample_folder + prefix + "{sample}" + readreverse_postfix,
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
        """--- Trimming ---"""
    shell:
        """trimmomatic PE {input.forward_strand} {input.reverse_strand} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired}   ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"""


#rule multiqc:
#    input: 
#        fastqc_done = expand(sample_work_path + "/fastqc/{sample}_fastqc.done", sample=MERGED_SAMPLES),
#        fastq_screen_done = expand(sample_work_path + "/fastqscreen/{sample}" + "_fastqscreen_done.txt", sample=MERGED_SAMPLES),
#    output:
#        qc_done = sample_work_path + "/multiqc/done.txt",
#    shell:
#        """ touch {output.qc_done} """
#
rule fastqc:
    input:
        forward_strand = raw_sample_folder + prefix + "{sample}" + readforward_postfix,
        reverse_strand = raw_sample_folder + prefix + "{sample}" + readreverse_postfix,
    output:
        fastqc_done = sample_work_path + "/fastqc/{sample}_fastqc_done.txt",
        #reverse = sample_work_path + "/fastqc/{sample}_2_fastqc.zip",
    params: 
        outdir= sample_work_path + "/fastqc/{sample}",
    conda:
        "../envs/fastqc.yaml"
    threads: 4
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """    
        fastqc --outdir  {params.outdir} --extract  -f fastq {input.forward_strand} {input.reverse_strand}
        touch {output.fastqc_done}
        """

#rule count reads
##This rule will count the number of reads in the fasta file and the number of mapped reads after alignment
##This will then be put into the multiqc quality repor that I am making
#The reason we are doing this as opposed to just looking at the percentage that mapped is because sometimes there arre reads that trimmomatic tosses out 
#because these specific reads are not able to be paired with each other and so trimmomatic throws them away

rule fastqscreen:
    input:
        forward_strand = raw_sample_folder + prefix + "{sample}" + readforward_postfix,
        reverse_strand = raw_sample_folder + prefix + "{sample}" + readreverse_postfix,
    output:
        #out_forward = sample_work_path + "/fastqscreen/{sample}/" + prefix + "{sample}" + readforward_postfix + "_screen.html",
        #out_reverse = sample_work_path + "/fastqscreen/{sample}/" + prefix + "{sample}" + readreverse_postfix + "_screen.html",
        fastq_screen_done = sample_work_path + "/fastqscreen/{sample}" + "_fastqscreen_done.txt",

    params:
        conf = fastq_screen_conf, 
        out_dir = sample_work_path + "/fastqscreen/{sample}"
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """
        fastq_screen --aligner bowtie2 --conf {params.conf} --outdir {params.out_dir} {input.forward_strand}
        fastq_screen --aligner bowtie2 --conf {params.conf} --outdir {params.out_dir} {input.reverse_strand}
        touch {output.fastq_screen_done}
        """



