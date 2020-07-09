# Snakemake workflow: atac_seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Rowan Callahan (@callahro)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/atac_seq/releases).
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

Installing the workflow involves cloning this directory and ensuring that you have Snakemake added to your path. To ensure that you have Snakemake installed you can run `which snakemake` and it will tell you which snakemake you have installed. Once you have cloned the workflow and checked your path you need to ensure that you have correctly formatted your configuration file. 



#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.
alternatively you can create your own configuration file and specify it when running snakemake with `snakemake --configfile my_config_file.yaml`
viewing the example configuration file will specify what values are needed where.


#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

This pipeline requires the usage of conda to run as all of its external dependencies are installed using conda.
its important that every time the pipeline is run that `--use-conda` is also called.

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/atac_seq) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


