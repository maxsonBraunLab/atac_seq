# atac_seq

[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)
![ci/cd status](https://github.com/maxsonBraunLab/atac_seq/actions/workflows/test.yaml/badge.svg)
[![snakemake minimum](https://img.shields.io/badge/snakemake->=5.32-<COLOR>.svg)](https://shields.io/)
![Maintainer](https://img.shields.io/badge/maintainer-gartician-blue)
![Maintainer](https://img.shields.io/badge/maintainer-bioThai-blue)

Process and analyze your PE ATAC-Seq datasets

## 1. Prepare your work environment

```bash
# clone this repo to a new working directory
git clone https://github.com/maxsonBraunLab/atac_seq.git

# cd into atac_seq-2.0 and make new dir for your FASTQ files
mkdir -p data/raw

# cd into 'data/raw' and symlink your FASTQ files
ln -s /absolute/path/to/files/condition1_R1.fastq.gz .
ln -s /absolute/path/to/files/condition1_R2.fastq.gz .

# rename symlinks to match {condition}_{replicate}_{read}.fastq.gz
mv condition1_R2 condition_1_R1.fastq.gz
mv condition1_R2 condition_1_R2.fastq.gz

# make scripts executable
# run chmod from main pipeline directory (where Snakefile is)
chmod +x scripts/*.py scripts/*.sh *.sh
```

Double check the symlinks match the required input format: `{condition}_{replicate}_{dir}.fastq.gz`.

## 2. Prepare your conda environment

Continue forward if you don't have a conda environment with snakemake installed

```bash
# while using base conda env, create your snakemake environment
conda install -c conda-forge mamba # installs into base env
mamba create -c conda-forge -c bioconda -n snakemake snakemake # installs snakemake into new env

conda activate snakemake
```

Make sure to also install plotly in the snakemake environment with `conda install -c plotly plotly` 

## 3. Prepare your pipeline configuration

Edit the `config/config.yaml` file to specify which organism to use and other pipeline parameters.

Run the `scripts/make_differential_configs.py` script from the main pipeline directory (where the Snakefile is located) to generate a `config/deseq2_config.tsv` and a `config/diffbind_config.csv` file:

```bash
# if snakemake env is not yet activated, then activate it
conda activate snakemake

# run python script
python scripts/make_differential_configs.py
```

The `config/deseq2_config.tsv` file specifies which replicates belong with which condition in DESeq2. Each column must be separated by a tab. The `config/diffbind_config.csv` file specifies sample names, sample conditions, and BAM file paths to use for Diffbind. 

Double-check that the `deseq2_config.tsv` and `diffbind_config.csv` files are properly formatted. Corresponding example files are provided in the `config` folder.

## 4. Set up SLURM integration (for batch jobs)

Do this step if are running the pipeline as a batch job and don't yet have a [SLURM profile](https://github.com/Snakemake-Profiles/slurm) set up.

Download the `slurm` folder from the maxsonBraunLab [repository](https://github.com/maxsonBraunLab/slurm) and copy the entire thing to `~/.config/snakemake`. 

Your file configuration for SLURM should be as follows:
```
~/.config/snakemake/slurm/<files>
```

Change the file permissions for the scripts in the `slurm` folder so that they are executable. To do this, run:
```
chmod +x ~/.config/snakemake/slurm/slurm*
```

## 5. Run the pipeline

You can run the pipeline using an interactive node like this:

```bash
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j 20 --use-conda
```

This is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process.

To run the pipeline using batch mode use the following command:

```bash
sbatch run_pipeline_conda.sh
```

Additional setup instructions are provided in the `run_pipeline_conda.sh` script. Make sure to follow these setup instructions prior to running the pipeline. 

For users running the pipeline in batch mode, `run_pipeline_conda.sh` is a wrapper script that contains the following command:

```bash
snakemake -j $num_jobs --verbose --use-conda --conda-prefix $CONDA_PREFIX_1/env --cluster-config cluster.yaml --profile slurm
```

This will submit up to `$num_jobs` jobs (This number can be changed in the wrapper script) to ARC servers and is appropriate for running computationally-intensive programs (read aligning, peak calling, finding consensus peaks, calculating differentially open chromatin regions).


## Pipeline Summary

### Assumptions

* When making the BWA index, chromosome names must be prefixed with 'chr' like 'chr1', 'chr2', 'chrX', and 'chrM'.

### Inputs

* Reads in this specific format: **{condition}\_{replicate}\_{dir}.fastq.gz**
    * `Condition` = experimental treatment such as: 'MOLM24D', 'SETBP1_CSF3R_Mutant' , 'SETBP1_CSF3R_Control_8hr'. Multiple underscores, text / number mixing is OK. 
    * `Replicate` = biological replicate. acceptable values are integers >= 1.
    * `Dir` = read direction. Acceptable values are ['R1', 'R2'].
    * Reads must be placed in the `data/raw` directory.
* Properly configured `config.yaml` file

### Outputs

All of the following are in the `data` directory.

* MultiQC report `data/multiqc/multiqc_report.html` that summarizes the following results:
    * sequencing quality from fastqc
    * sequencing quality and read trimming from fastp
    * alignment quality from bowtie2
    * library contamination from fastq_screen
    * library complexity from preseq

* FRiP and fragment length distribution `data/frip.html` `data/fraglen.html`. FRiP is more like Fraction of Reads in Consensus Peaks.

* Consensus peaks among _n_ replicates per condition where _n_ is configurable. A peak is considered a consensus peak when it is in >= _n_ number of samples per condition. 

    * For example, peak1 appears in 2 out of 3 replicates in condition1. If n = 2, then peak1 is considered a consensus peak, even if it does not appear in other conditions. 
    * For example, peak2 appears in 1 out of 3 replicates in all conditions. If n = 2, then peak2 is not a consensus peak. 

* Raw counts table of peaks (rows = intervals, cols = samples) `data/counts/counts_table.txt`

* Bigwig files using CPM normalization `data/bigwig`

* Differentially open chromatin regions for all unique combinations of conditions. Instead of specifying contrasts explicitly, the pipeline will assess all unique combinations of conditions.

    * PCA plot of all samples `data/deseq2/sample_PCA.png`
    * DESeq2-normalized and ln(DESeq2-normalized) counts `data/deseq2/norm_counts.txt`, `data/deseq2/log_norm_counts.txt`
    * Significant peaks split by up and down-regulation `data/deseq2/{contrast}/{contrast}-[up_sig|down_sig].bed`
    * Heatmap of top 50 most differential regions `data/deseq2/{contrast}/{contrast}-heatmap.pdf` 

* HOMER analysis of up and down differential peaks per contrast `data/homer/{contrast}-{up|down}`

    * HOMER will only run if there is >= 10 differential up/down peaks

## Methods

![](rulegraph.svg).



