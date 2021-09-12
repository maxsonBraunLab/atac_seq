# readme

`.test` is the testing environment for `atac_seq`. It is used with Continuous Integration in mind. Below are the instructions to prepare the data for CI/CD.

## 0. Download Data and Set Up Testing Environment

### Download data to `data/raw`

```bash
mkdir -p data/raw
cd data/raw

# the 1 and 2 really mean different lanes, not bio replicate. but doesn't matter in test case.

# young adrenal glands
wget -q https://www.encodeproject.org/files/ENCFF088ONP/@@download/ENCFF088ONP.fastq.gz &
wget -q https://www.encodeproject.org/files/ENCFF011ANA/@@download/ENCFF011ANA.fastq.gz &

wget -q https://www.encodeproject.org/files/ENCFF097DWT/@@download/ENCFF097DWT.fastq.gz &
wget -q https://www.encodeproject.org/files/ENCFF501VCA/@@download/ENCFF501VCA.fastq.gz &

# old adrenal glands
wget -q https://www.encodeproject.org/files/ENCFF106DBO/@@download/ENCFF106DBO.fastq.gz &
wget -q https://www.encodeproject.org/files/ENCFF093HFA/@@download/ENCFF093HFA.fastq.gz &

wget -q https://www.encodeproject.org/files/ENCFF126VJY/@@download/ENCFF126VJY.fastq.gz &
wget -q https://www.encodeproject.org/files/ENCFF817MYO/@@download/ENCFF817MYO.fastq.gz &

# subset and rename
zcat ENCFF088ONP.fastq.gz | head -4000000 | gzip > AdrenalGland16_1_R1.fastq.gz &
zcat ENCFF011ANA.fastq.gz | head -4000000 | gzip > AdrenalGland16_1_R2.fastq.gz &
zcat ENCFF097DWT.fastq.gz | head -4000000 | gzip > AdrenalGland16_2_R1.fastq.gz &
zcat ENCFF501VCA.fastq.gz | head -4000000 | gzip > AdrenalGland16_2_R2.fastq.gz &

zcat ENCFF106DBO.fastq.gz | head -4000000 | gzip > AdrenalGland59_1_R1.fastq.gz &
zcat ENCFF093HFA.fastq.gz | head -4000000 | gzip > AdrenalGland59_1_R2.fastq.gz &
zcat ENCFF126VJY.fastq.gz | head -4000000 | gzip > AdrenalGland59_2_R1.fastq.gz &
zcat ENCFF817MYO.fastq.gz | head -4000000 | gzip > AdrenalGland59_2_R2.fastq.gz &

mv ENCFF088ONP.fastq.gz ENCFF088ONP.fastq.gz.norun
mv ENCFF011ANA.fastq.gz ENCFF011ANA.fastq.gz.norun
mv ENCFF097DWT.fastq.gz ENCFF097DWT.fastq.gz.norun
mv ENCFF501VCA.fastq.gz ENCFF501VCA.fastq.gz.norun

mv ENCFF106DBO.fastq.gz ENCFF106DBO.fastq.gz.norun
mv ENCFF093HFA.fastq.gz ENCFF093HFA.fastq.gz.norun
mv ENCFF126VJY.fastq.gz ENCFF126VJY.fastq.gz.norun
mv ENCFF817MYO.fastq.gz ENCFF817MYO.fastq.gz.norun

```

### Download and index Chr1 of hg38

```bash
mkdir testdata

wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.1.fa testdata
bwa index -p testdata/hg38_chr1 Homo_sapiens.GRCh38.dna.chromosome.1.fa

```

### Fill out the `config.yaml`

Put genomes and banlist to `testdata`

```yaml
# pipeline config ---------------------------------------------------------------------------------

# please adjust content of these files, though their paths should be constants.
FASTQ_SCREEN_CONFIG: "data/config/fastq_screen_config.tsv"
MULTIQC_CONFIG: "data/config/multiqc_config.yaml"
DESEQ2_CONFIG: "data/config/deseq2_config.tsv"

# read alignment + preprocessing ------------------------------------------------------------------

GENOME: "testdata/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
ASSEMBLY: "hg38" # [hg38, mm10]

# Forbidden regions of the genome for read alignment and analysis
BANLIST: "testdata/hg38.blacklist.v2.bed"

# peak calling ------------------------------------------------------------------------------------

# number of times a peak needs to appear across samples in one condition. helps build consensus peak.
N_INTERSECTS: 2

# size of the genome for macs2 peak calling
GSIZE: "2.7e9"

# differential analysis ---------------------------------------------------------------------------
# filter for results with pval less than <padj_cutoff>
padj_cutoff: 0.05

```

### Miscellaneous

Download _saccharomyces cerevesiae_ bowtie2 index for fastq_screen into `testdata/screen`.

## 2. Run pipeline

From this point on, just follow the original `readme.md` file to process the reads. This includes filling out the `config.yaml`. From here, continuous integration will load up Conda packages and execute the pipeline.

## References

S. Triana, D. Vonficht, et al. [Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states](https://www.biorxiv.org/content/10.1101/2021.03.18.435922v1.full#disqus_thread).