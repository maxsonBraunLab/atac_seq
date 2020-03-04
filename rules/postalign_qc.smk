rule makeBigwig:
    input:
        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
        deduplicated = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
    output:
        bigwig = sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_smooth.bw",
        bigwig_rough = sample_work_path + "/fully_filtered/{merged_sample}_tracks_rpkm_5window_rough.bw",

    conda:
        "../envs/deeptools.yaml"
    threads: 8 
    message:
        """--- making bigwig---"""
    shell:
        """
        bamCoverage --bam {input.quality} -o {output.bigwig} --numberOfProcessors {threads} --skipNAs --binSize 5 --smoothLength 15  --normalizeUsing RPKM
        bamCoverage --bam {input.quality} -o {output.bigwig_rough} --numberOfProcessors {threads} --skipNAs --binSize 5 --normalizeUsing RPKM

        """

rule shiftReads:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality_shiftedReads.bam"
    conda:
        "../envs/sambamba.yaml"
    threads: 2
    message:
        """--- shifting reads---"""
    script:
        "../scripts/shift_reads.R" 

rule index:
    input:
        deduplicated = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam",
        full_bam = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
    output:
        indexed = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam.bai",
        indexed_full = sample_work_path + "/bamfiles/{merged_sample}_merged.bam.bai",
    conda:
        "../envs/sambamba.yaml"
    threads: 4 
    message:
        """--- indexing reads---"""
    shell:
        """
        sambamba index -t {threads} {input.full_bam} {output.indexed_full}
        sambamba index -t {threads} {input.deduplicated} {output.indexed}
        
        """

rule countReads:
     input:
        merged = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
        dedup = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
        no_mito = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam",
        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam"
     output:
        readsfile = sample_work_path + "/fully_filtered/{merged_sample}_read_depths.csv"
     conda:
        "../envs/samtools_env.yaml"
     threads: 4
     message:
        """--- counting reads ---"""
     shell:
        """
        samtools view -@ 4 -c {input.merged}| xargs -I{{}} echo "{wildcards.merged_sample},starting_reads,"{{}}$'\n' >> {output.readsfile}
        samtools view -@ 4 -c {input.no_mito}| xargs -I{{}} echo "{wildcards.merged_sample},no_ChrM,"{{}}$'\n' >> {output.readsfile}
        samtools view -@ 4 -c {input.dedup}| xargs -I{{}} echo "{wildcards.merged_sample},no_duplicates,"{{}}$'\n'  >> {output.readsfile}
        samtools view -@ 4 -c {input.quality}| xargs -I{{}} echo "{wildcards.merged_sample},satisfy_quality,"{{}}$'\n'  >> {output.readsfile}
        """


rule samtoolsQuality:
     input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam",
     output:
        quality = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup_quality.bam"
     conda:
        "../envs/samtools_env.yaml"
     threads: 4
     message:
        """--- deduplicating reads---"""
     shell:
        """samtools view -h -@ 4 -b -F 1804 {input.bamfile} > {output.quality}"""
        # Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
        # Retain properly paired reads -f 2
 
rule dedup:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
    output:
        deduplicated = sample_work_path + "/bamfiles/{merged_sample}_rmChrM_dedup.bam"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
        """--- deduplicating reads---"""
    shell:
         """sambamba markdup -t {threads} {input.bamfile} {output.deduplicated}"""

rule removeMitochondrial:
    input:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_sorted.bam", #change this to _merged if merging lanes
    output:
        bamfile = sample_work_path + "/bamfiles/{merged_sample}_rmChrM.bam",
    conda:
        "../envs/samtools_env.yaml"
    message:
        """--- Removing mitochondrial reads---"""
    log:
        totalreads = sample_work_path + "/bamfiles/{merged_sample}_readswmito.txt"
    threads: 4 
    shell:
        """
        samtools view -@ {threads} -c {input.bamfile} > {log.totalreads}
        samtools view -h -@ {threads} {input.bamfile}| grep -v chrM| samtools sort -@ {threads} -O BAM > {output.bamfile}
        """

rule mergeBam:
    input:
        lane1 = sample_work_path + "/bamfiles/{merged_sample}" + LANES[0] + "_sorted.bam",
        lane2 = sample_work_path + "/bamfiles/{merged_sample}" + LANES[1] + "_sorted.bam",
    output:
        merged_bamfile = sample_work_path + "/bamfiles/{merged_sample}_merged.bam",
    conda:
        "../envs/samtools_env.yaml"
    threads: 8 
    message:
         """--- sorting the bamfile ---"""
    shell:
         """samtools merge -@ {threads} {output.merged_bamfile} {input.lane1} {input.lane2}"""



