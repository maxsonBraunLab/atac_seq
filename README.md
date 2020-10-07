# atac_seq2

Process and analyze your ATAC-Seq datasets

# 1. Prepare your work environment

```bash
# clone this repo to a new working directory
git clone git@github.com:maxsonBraunLab/atac_seq-2.0.git

# cd into atac_seq-2.0 and make new dir for your FASTQ files
mkdir -p samples/raw

# cd into 'samples/raw' and symlink your FASTQ files
ln -s /absolute/path/to/files/condition1_replicate1_R1.fastq.gz .
ln -s /absolute/path/to/files/condition1_replicate1_R2.fastq.gz .
```



# 2. Prepare your conda environment

Continue forward if you don't have a conda environment with snakemake installed

```bash
# while using base conda env, create your snakemake environment
conda install -c conda-forge mamba # installs into base env
mamba create -c conda-forge -c bioconda -n snakemake snakemake # installs snakemake into new env

conda activate snakemake
```

# 3. Run the pipeline

You can run the pipeline using an interactive node like this:

```bash
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j 20 --use-conda
```

This is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process. 



You can run the pipeline via batch mode like this:

```bash
snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -J {cluster.J} -c {cluster.c} --mem={cluster.mem}" -s Snakefile
```

This will submit up to 64 jobs to exacloud servers and is appropriate for running computationally-intensive programs (read aligning, peak calling, calculating differentially open chromatin regions).



# Pipeline Summary



## Inputs

* Reads in this specific format: **{condition}\_{replicate}\_{dir}.fastq.gz**
  * `Condition` = experimental treatment such as: 'MOLM24D', 'SETBP1_CSF3R_Mutant' , 'SETBP1_CSF3R_Control_8hr'. Multiple underscores, text / number mixing is OK. 
  * `Replicate` = biological replicate. acceptable values are integers >= 1.
  * `Dir` = read direction. Acceptable values are ['R1', 'R2'].
  * Reads must be placed in the `samples/raw` directory.
* Adapter file in FASTA format for adapter trimming.
* Reference genome in FASTA format.

## Outputs

* Quality Control
  * Fragment length distribution plot
  * Fraction of Reads in Peaks (FRiP) per sample 
  * PCA of all replicates
* Quality Table of number of reads after removing mito, duplicates, poorly mapping
* Counts table of peaks
  * Pre and post-downsampling to sample with the lowest coverage 
* Fraction of Reads in Peaks (FRIP) per sample
* Consensus peaks among _n_ replicates (_n_ is configurable) 
* Read pileup tracks in bigwig format
* Differentially open chromatin regions
* Commonly used data (tracks, QC metrics, counts table) are in `data` directory.

## Methods

![](rulegraph.svg).

# References


