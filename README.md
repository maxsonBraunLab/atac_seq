# atac_seq2

Process and analyze your PE ATAC-Seq datasets

# 1. Prepare your work environment

```bash
# clone this repo to a new working directory
git clone git@github.com:maxsonBraunLab/atac_seq2.git

# cd into atac_seq-2.0 and make new dir for your FASTQ files
mkdir -p samples/raw

# cd into 'samples/raw' and symlink your FASTQ files
ln -s /absolute/path/to/files/condition1_replicate1_R1.fastq.gz .
ln -s /absolute/path/to/files/condition1_replicate1_R2.fastq.gz .
```

After the files are symlinked, rename the symlinks to match the required input format: `{condition}_{replicate}_{dir}.fastq.gz` (more information in pipeline summary below).

# 2. Prepare your conda environment

Continue forward if you don't have a conda environment with snakemake installed

```bash
# while using base conda env, create your snakemake environment
conda install -c conda-forge mamba # installs into base env
mamba create -c conda-forge -c bioconda -n snakemake snakemake # installs snakemake into new env

conda activate snakemake
```

Make sure to also install plotly in the snakemake environment with `conda install -c plotly plotly` 

# 3. Prepare your pipeline configuration

Edit the `config.yaml` file to specify which organism to use and other pipeline parameters.

Edit the `config/metadata.csv` file to specify which replicates belong with which condition in DESeq2.

# 4. Run the pipeline

Copy the entire slurm folder into your home directory `~/.config/snakemake` and then delete the local copy.

You can run the pipeline using an interactive node like this:

```bash
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j 20 --use-conda
```

This is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process. 



You can run the pipeline via batch mode like this:

```bash
snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --cluster-config cluster.yaml --profile slurm

# if the above does not work, try the command below
snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -J {cluster.J} -c {cluster.c} --mem={cluster.mem}" -s Snakefile

```

This will submit up to 64 jobs to exacloud servers and is appropriate for running computationally-intensive programs (read aligning, peak calling, finding consensus peaks, calculating differentially open chromatin regions).

# Pipeline Summary



## Inputs

* Reads in this specific format: **{condition}\_{replicate}\_{dir}.fastq.gz**
  * `Condition` = experimental treatment such as: 'MOLM24D', 'SETBP1_CSF3R_Mutant' , 'SETBP1_CSF3R_Control_8hr'. Multiple underscores, text / number mixing is OK. 
  * `Replicate` = biological replicate. acceptable values are integers >= 1.
  * `Dir` = read direction. Acceptable values are ['R1', 'R2'].
  * Reads must be placed in the `samples/raw` directory.
* Adapter file in FASTA format for adapter trimming. This is provided.
* Reference genome in FASTA format.

## Outputs

All of the following are in the `data` directory.

* Quality Control plots in HTML format

  * Fragment length distribution plot `fragment_len_dist.html`
  * Fraction of Reads in Peaks (FRiP) per sample `frip.html`
  * Report of the data with important metrics at alignment, peak calling, and DE steps `essential_report.html`. This includes 

* Raw counts table of peaks (rows are intervals, columns are samples) `counts_table.txt`

* Read pileup tracks in bigwig format using CPM normalization `data/bw/*.bw`

* Differentially open chromatin regions for all unique combinations of conditions. Instead of specifying contrasts explicitly, the pipeline will assess all unique combinations of conditions.

  * PCA plot of all samples `data/de/sample_PCA.png`
  * DESeq2-normalized and ln(DESeq2-normalized) counts `de/norm_counts.txt`, `de/log_norm_counts.txt`
  * Significant peaks split by up and down-regulation `de/{contrast}/{contrast}-[up_sig|down_sig].bed`
  * Heatmap of top 50 most differential regions `de/{contrast}/{contrast}-heatmap.pdf` 

* Consensus peaks among _n_ replicates per condition (_n_ is configurable). This is in the `samples` folder. A peak is considered a consensus peak when it is in >= _n_ number of samples per condition. 

  * For example, peak1 appears in 2 out of 3 replicates in condition1. If n = 2, then peak1 is considered a consensus peak, even if it does not appear in other conditions. 
  * For example, peak2 appears in 1 out of 3 replicates in all conditions. If n = 2, then peak2 is not a consensus peak. 

  # Methods

![](rulegraph.svg).

# References
