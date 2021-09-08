# atac_seq

Process and analyze your PE ATAC-Seq datasets

## 1. Prepare your work environment

```bash
# clone this repo to a new working directory
git clone git@github.com:maxsonBraunLab/atac_seq.git

# cd into atac_seq-2.0 and make new dir for your FASTQ files
mkdir -p data/raw

# cd into 'data/raw' and symlink your FASTQ files
ln -s /absolute/path/to/files/sample_R1.fastq.gz .
ln -s /absolute/path/to/files/sample_R2.fastq.gz .

# rename symlinks to match {condition}_{replicate}_{read}.fastq.gz
mv sample_R1.fastq.gz condition_replicate_R1.fastq.gz
mv sample_R2.fastq.gz condition_replicate_R2.fastq.gz
```

Double check that the symlinks match the required input format: `{condition}_{replicate}_{read}.fastq.gz`. Condition can contain text, numbers, and dashes; replicate must contain integers; read is either R1 or R2.

## 2. Prepare your conda environment

Continue forward if you don't have a conda environment with a clean install of snakemake.

```bash
# while using base conda env, create your snakemake environment
conda install -c conda-forge mamba # installs mamba into base env
mamba create -c conda-forge -c bioconda -n snakemake snakemake # installs snakemake into new env with mamba
conda activate snakemake
```

Make sure to also install plotly in the snakemake environment with `conda install -c plotly plotly` 

## 3. Prepare your pipeline configuration

Edit the `config.yaml` file to specify which organism to use and other pipeline parameters.

Run the following command to generate differential config files for DESeq2 and DiffBind:

```bash
$ python scripts/make_differential_configs.py
```

Then double check the `data/config/deseq2_config.tsv|diffbind_config.csv` files. You can rename the condition columns to change the output in DESeq2/DiffBind. 

## 4. Run the pipeline

To configure SLURM snake make profile with OHSU's cluster, copy the entire `slurm` folder into your home directory `~/.config/snakemake` and then delete the local copy.

You can run the pipeline using batch mode like this:

```bash
snakemake -j 64 --use-conda --cluster-config cluster.yaml --profile slurm --restart-times 3
```

This will submit up to 64 jobs to exacloud servers and is appropriate for running computationally-intensive programs (read aligning, peak calling, finding consensus peaks, calculating differentially open chromatin regions).

---

You can also run the pipeline interactively like this:

```bash
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j 20 --use-conda --restart-times 3
```

This is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process.

## Output

All of the following are in the `data` directory.

```
data/banlist ------- post-processed bam files
data/bigwig -------- tracks for IGV
data/config -------- config file for multiqc, DE testing, fastq_screen
data/counts -------- consensus peaks + counts tables per sample
data/deseq2 -------- results of DESeq2
data/diffbind ------ results of DiffBind
data/fastp --------- trimmed and QC'd reads
data/fastqc -------- results of FASTQC
data/fastq_screen -- results of fastqc_screen
data/logs ---------- log files
data/macs2 --------- results of MACS2
data/multicov ------ counts table containing all samples
data/multiqc ------- quality control report
data/preseq -------- results of Preseq
data/raw ----------- raw data
```

## Pipeline Details

FRiP is actually Fraction of Reads in Consensus Peaks.

Consensus peaks is the union of peaks that appear in _n_ replicates in at least one condition where _n_ is configurable. 

* For example, peak1 appears in 2 out of 3 replicates in condition1. If n = 2, then peak1 is considered a consensus peak, even if it does not appear in other conditions. 
* For example, peak2 appears in 1 out of 3 replicates in all conditions. If n = 2, then peak2 is not a consensus peak. 

Bigwig files are made with CPM normalization to compare between samples without sequencing depth bias.

DESeq2 and DiffBind will test for differential open chromatin regions using every unique combinations of conditions. No need to specify your contrasts anymore, we took care of that! DESeq2 outputs contain:

* PCA plot of all samples `data/deseq2/sample_PCA.png`

* DESeq2-normalized and ln(DESeq2-normalized) counts `deseq2/norm_counts.txt`, `deseq2/log_norm_counts.txt`
* Significant peaks split by up and down-regulation `deseq2/{contrast}/{contrast}-[up_sig|down_sig].bed`
* Heatmap of top 50 most differential regions `deseq2/{contrast}/{contrast}-heatmap.pdf` 

DESeq2 normalizes data by the reads in peaks and assumes that most features (intervals in this case) do not change. DiffBind normalizes data by sequencing depth. We find that DESeq2 works well for finding differential peaks, while DiffBind works better to identify global changes in open chromatin regions. While the biological patterns will differ based on experiment, we want to give both options for the user to explore the nuance in their results..

HOMER motif analysis of up and down differential peaks will only be run if there is > 10 peaks. This rule will submit each HOMER run to SLURM, so make sure you are using a SLURM-based HPC. It's also easy to configure it to something else if needed. It currently only supports the output of DESeq2. 

## Methods

![](rulegraph.svg).
