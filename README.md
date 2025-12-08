# MetaDetector

## Summary

MetaDetector is tool designed for identification and classification of metagenomic sequencing data. \
This tool is primarily written with shell scripting and implemented for high performance computing clusters using SLURM workload manager. This README will provide the user with instructions for installation of proper tools and databases to run MetaDetector.


## Files

 - The bin directory contains the scripts required to run the pipeline
 - This README and databases/README should be followed prior to use of MetaDetector
 - The md.yml file is used to create the conda environment

## SLURM
 - MetaDetector is designed to be used with SLURM Workload Manager
 - SLURM installation instructions can be found here: https://slurm.schedmd.com/quickstart_admin.html 

## Installation
Setting up the environment: \
    `conda create -n md --file /path/to/md.yml` \
    `conda activate md`
 - Creating and activating this environment should provide all of the appropriate packages required for setting up and running Metadetector
 - If not already installed, Mamba can be installed here: https://conda-forge.org/download/

Installing Databases:
 - Please refer to databases/README for database installation
 - The previous step, setting up the conda environment, must be completed prior to database installation
 - The steps outlined in the databases/README should be followed and completed prior to completing any further steps in the README

Setting Paths:
 - Paths to various directories used in MetaDetector are saved in the `md.config` file
 - The md.config file should be located in the same directory as the scripts used to run MetaDetector
 - Generic paths are currently in the file, most of which to a "databases" directory and subdirectories, but you are free to organize this as you see fit
 - **NOTE**: Other scripts (helper, update, and database) have hard-coded paths that will need to be manually modified

## Running MetaDetector
#### Create a readlist:
 - Create a text file with paths to fastq files. For example: \
 `ls -d -w 1 /path/to/fastq/files/*fastq.gz > /full/path/to/readlist`
 - Ex. A samplesheet containing reads from both illumina and minion sequencing platforms
 ```
 /path/to/Sample_A_R1.fastq.gz
 /path/to/Sample_A_R2.fastq.gz
 /path/to/Sample_A.fastq.gz
 ```
 - Illumina reads require paths to both R1 and R2 files, otherwise the pipeline will not run
 - If running in "hybrid" mode the sample names must contain the same basename (i.e. "Sample_A" in the example above)

#### Pipeline Options: 
 `[-a assembly_flag] - (default "m") "a" for all assembly stages, "i" for isolate assembly, "m" for metagenomic assembly, "t" for isolate with targeted reference matches, "v" for metaviral  assembly` \
 `[-b bbmap options] - bbmap options to run for read mapping for host and contaminant removal` \
 `[-d host_db for host removal] - /path/to/host/fasta file for removal or targeted reference matching` \
 `[-h help]` \
 `[-j assembly_type] - (default "s") type of assembly according to reads e.g. "s" (short), "h" (hybrid), "l" (long)` \
 `[-m memory in GB] - amount of memory to allocate on slurm partition (in GB)` \
 `[-o output_directory] - name of output directory WITHOUT the full path` \
 `[-r readlist] - file containing paths to fastq files` \
 `[-t threads] - number of threads to allocate on SLURM and multithreaded processes` \

#### Run the pipeline:
 - Run MetaDetector on a set of illumina sequencing data:  
 `bash queue.sh -r /full/path/to/readlist -o /full/path/to/output_dir `
 - Run MetaDetector on a set of Minion sequencing data:  
 `bash queue.sh -r /full/path/to/readlist -o /full/path/to/output_dir -j "l"`
 - Run MetaDetector on a mix of Illumina and Minion sequencing data and run in isolate mode:  
 `bash queue.sh -r /full/path/to/readlist -o /full/path/to/output_dir -j "h" -a "i"`

## Outputs
 - Primary outputs from MetaDetector are `[SampleName]_contigs.daa` and `[SampleName]_reads.daa` files. These files can be opened with Megan7 and further analysis may be performed there. 
 - Additional outputs include:
   - Trimmed reads
   - Host and contaminant removed reads
   - Assembly outputs - contigs and assembly graphs
   - Alignment of reads against contigs
   - Quast report of assembly
 - Supplemental files:
   - Scripts for each stage of the pipeline
   - Log output and error files 
   - Status log files to bookmark successfully completed stages

## Disclaimer
This work was supported by work unit number A1714.

The views expressed in this article reflect the results of research conducted by the author and do not necessarily reflect the official policy or position of the Department of the Navy, Department of Defense, nor the U.S Government. Several authors are military service members or federal employees of the United States government.

This work was prepared as part of their official duties. Title 17 U.S.C. 105 provides that ‘copyright protection under this title is not available for any work of the United States Government.’ Title 17 U.S.C. 101 defines a U.S. Government work as work prepared by a military service member or employee of the U.S. Government as part of that person's official duties.
