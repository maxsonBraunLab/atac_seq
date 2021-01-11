rule chip_screen:
    input:
        motifs = "samples/macs/reads_catalog_intervals_nodownsample.bed",
        chip = "/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/tf_chip.bed.gz"
    output:
        "data/chip_screen/{db}/{cell_line}/consensus_peaks_tf_chip.bed.gz"
    conda:
        "../envs/frip.yaml"
    shell:
        "bedtools intersect -wa -wb -f 0.50 -a {input.motifs} -b {input.chip} | gzip > {output}"

rule chip_fisher:
    input:
        peaks = "samples/macs/reads_catalog_intervals_nodownsample.bed",
        chip = "/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/tf_chip.bed.gz"
    output:
        "data/chip_screen/{db}/{cell_line}/consensus_peaks_tf_chip.out"
    params:
        chrom_sizes = config["chrom_sizes"]
    conda:
        "../envs/frip.yaml"
    shell:
        "bedtools merge -i {input.chip} | bedtools fisher -f 0.50 -a {input.peaks} -b stdin -g {params.chrom_sizes} > {output}"
# merge all tf's in db into one interval to reduce multiple intersects + account for different size db's.

rule frip_plot:
    input:
        expand("samples/frip/{sample}_stats.txt", sample = SAMPLE)
    output:
        "data/frip.html"
    run:
        pd.options.plotting.backend = "plotly"
        # import data. row = total_reads, reads_in_peaks, and ratio. cols = samples. 
        df = pd.concat([ pd.read_csv(i) for i in sorted(input) ], axis = 1)
        df = df.rename(index={0: 'total_reads', 1: 'reads_in_peaks'})
        df.loc['ratio'] = df.loc['reads_in_peaks'] / df.loc['total_reads']
        # plot graph. plot ratio as bottom as percent, and plot to max value of 1.
        fig = go.Figure(data=[
            go.Bar(name='FRiP', x=df.columns, y=df.loc['ratio'], marker_color='rgb(255,201,57)'),
            go.Bar(name='fraction_reads', x=df.columns, y= ([1] * df.shape[1]) - df.loc['ratio'], marker_color='rgb(0,39,118)')])
        # Change the bar mode
        fig.update_layout(barmode='stack', title='Fraction of Reads in Peaks by Sample', xaxis_tickfont_size=14,
            yaxis=dict(title='Fraction of reads in peaks', titlefont_size=16, tickfont_size=14),
            xaxis=dict(title='Samples'))
        fig.write_html(str(output))

rule frip_count:
    input:
        bamfile = "samples/bamfiles/filtered/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam",
        peaks = "samples/macs/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak"
    output:
        "samples/frip/{sample}_stats.txt"
    conda:
        "../envs/frip.yaml"
    shell:
        """
        rip=$(bedtools sort -i {input.peaks} | bedtools merge -i stdin | bedtools intersect -u -a {input.bamfile} -b stdin -ubam | samtools view -c); \
        counts=$(samtools view -c {input.bamfile}); \
        echo -e "{wildcards.sample}\n$counts\n$rip" > {output}
        """
# source: yiwei niu

#This rule merges all of the catalog counts from the individual samples into one catalog for all of the samples
rule counts_table:
    input:
        read_count = expand("samples/macs/counts/{sample}_counts.bed", sample=SAMPLE),
    output:
        read_count = "data/counts_table.txt"
    run:
        import pandas as pd
        dataframes = []
        for file in input.read_count:
            dataframes.append(pd.read_csv(file, sep='\t'))
        merged = pd.concat(dataframes, axis=1)
        merged = merged.loc[:,~merged.columns.duplicated()]   
        # merged = merged.sort_index(axis=1)
        merged.to_csv(output.read_count, header=True, index=False, sep='\t')

# counts per sample at consensus intervals
rule sample_counts:
    input:
        catalog = "samples/macs/reads_catalog_intervals_nodownsample.bed",
        bam = "samples/bamfiles/filtered/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam"
    output:
        "samples/macs/counts/{sample}_counts.bed"
    params:
        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + "{sample}"
    conda:
        "../envs/peaks_catalog.yaml"
    shell:
        "bedtools multicov -bams {input.bam} -bed {input.catalog} > {output}; sed -i '1i {params.header}' {output}"

#This rule has two outputs, one is the overall peak catalog and another is the peak catalog for each sample condition
#This can later be input into the bedtools commands to allow us to assign peaks in the sample condition catalogs to the 
#corresponding peak in the overall peak catalog so that we can run bedtools multiintersect
rule peak_catalog_no_downsmpl:
    input:
        peaks = expand("samples/macs/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak", sample=SAMPLE),
    output:
        reads_catalog_bed = "samples/macs/reads_catalog_intervals_nodownsample.bed",
        consensus_stats = "samples/macs/consensus_stats.txt"
    params:
        merged_bed = "samples/macs/all_broad_peaks_nodownsample.bed",
        present_in_number = config["n_intersects"],
        blacklist = config["blacklist_file"],
        genome = config["genome_name"],
        metadata = config["metadata_file"],
        regex = "", #This assumes that the sample names will exactly match 
        #What is on the metadata file. This feature should be removed as
        #its implementation is lacking
        script_file = "./scripts/computing_peak_catalog_by_replicate.R"
    conda:
        "../envs/peaks_catalog.yaml"
    threads: 8
    shell:
        """
        cat {input.peaks} > {params.merged_bed}
        Rscript --vanilla {params.script_file} {params.merged_bed} {params.blacklist} {params.present_in_number} {output.reads_catalog_bed} {params.genome} {params.metadata} {params.regex} 
        """

#This command runs MACS on the non downsampled bam files
rule MACS_no_downsmpl:
    input:
        bamfile = "samples/bamfiles/filtered/{sample}_rmChrM_dedup_quality_shiftedReads_sorted.bam"
    output:
        summits = "samples/macs/{sample}_macsout_nodownsample/{sample}_macs_peaks.broadPeak"
    conda:
        "../envs/MACS.yaml"
    threads: 4
    params:
        genome = config["genome_size"],
        name = "{sample}_macs",
        directory = "samples/macs/{sample}_macsout_nodownsample"
    message:
        """--- running MACS ---"""
    shell:
        """
        macs2 callpeak --treatment {input.bamfile} --name {params.name} --format BAMPE --gsize {params.genome} --outdir {params.directory} --broad
        """


