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

An example configuration file is included that has information around what indices and databases to use.
Some of these configurations are paths to databases or whole genomes. Some of these databases have already been downloaded by me.
If you are using Human sample data then you can probably use some of the defaults in the pipeline that are included in the file titled `Snakefile`.
This file also includes the final target outputs that will be created by the pipeline and can be edited if you want to include or not include certain final outputs from the pipeline.

**Metadata File configuration**
The metadata file is required to have two columns one titled SampleID and the other titled Condition
This file also must be tab separated. When you run this pipeline in a dry run it will create a list of the sample names that it has found and print them out
The script works by matching these sample names exactly to the SampleID strings that are present in the metadata file.
One exception is if you use the `remove_filepath_regex:` option in the configuration file. This will perform a string replacement on the sample names that are listed
before trying to match them to the SampleID strings. This can be useful if a large part of the filepath does not contain sample identifiers and is instead repeated.
In this case you can simply write the part of the filepath that you care about in SampleID and cut off the rest with `remove_filepath_regex:`.

An example of this is as follows:

lets say I have samples in my sample list called `['sample1_long_filepath_to_ignore', 'sample2_long_filepath_to_ignore']`
I can set my `remove_filepath_regex:` in my config file to `remove_filepath_regex = '_long_filepath_to_ignore'` .
This would mean that now the SampleID list is looking for matches that look like this `['sample1','sample2']` these matches are used to create a dictionary that assigns 
sample names to the condition that they are present in.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

This pipeline requires the usage of conda to run as all of its external dependencies are installed using conda.
its important that every time the pipeline is run that `--use-conda` is also called.

It is important that you remember to specity to use a profile path using `--profile /path/to/profile/folder`.
For different versions of Snakemake sometimes this will lead to an error saying it did not find a profile configuration file.
If this is the case then it is necessary that you provide the path to the profile configuration file with `--profile /path/to/profile/config.yaml`.

It is possible to submit your snakemake command with the sbatch directive but it is also possible to just use a terminal multiplexer like `screen` or `tmux` as the snakemake scheduler process does not take up that much cpu or memory. However, its important that your profile is set up so that it sends jobs to the queue and does not run these jobs on the head node.


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


