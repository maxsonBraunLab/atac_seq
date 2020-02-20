#!/bin/bash
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=10:00:00
#SBATCH -c 20
####  This command works so far so this should be exactly the same as what this sbatch script does
###   sbatch -c 12 --mem-per-cpu=2gb snakemake --use-conda --rerun-incomplete --profile /home/groups/MaxsonLab/callahro/atac_seq/profile/


snakemake --rerun-incomplete --use-conda --profile ./profile/

