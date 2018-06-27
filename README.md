# ATAC2GRN Project README

This is the Github for ATAC2GRN, a validated pipeline for running chromatin profiling assays. Here you'll find code that was used to produce the figures in the paper as well as pipelines and tools for validating your own pipelines.

## How to Get Started

If you're looking to take bits and pieces of this code for your own use, skip forward to that folder's discussion. Otherwise, here are two options to set up your machine with the appropriate prerequisites.

### Singularity File

Attached is the file

```
ATAC2GRN_Singularity
```

This is a file that specifies commands in an order to a program that will build a Linux system (Ubuntu 16.04, specifically). If you're comfortable using Singularity, create a container with this file on a system using the Ubuntu 16.04 OS (approximate size ~35GB). If you'd like to learn how to build a Singularity container using this file, a quick tutorial can be found [here](https://singularity.lbl.gov/quickstart).

We'd like to thank the NIH HPC for the inspiration they provided with these Singularity files.

### Using VirtualBox and a VDI File

If you'd like to build a machine in VirtualBox to run these pipelines, download the .vdi file from [our Google Drive](https://drive.google.com/open?id=1j-sO0CjyK-u95Y2ZPQIDz5zqH2KnvAnd).

Unpack the file with

```
tar -xvzf ATAC2GRN_VDI.tar.gz
```

Inside the folder will be a .vdi file, and this can be loaded into VirtualBox directly as an Ubuntu 16.04 system with all of the dependencies for this pipelines installed. We'd recommend allocating sizable hard drive space. We used 400 GB, though as a rule of thumb you should have at least 250% as much space as it takes to accomodate your raw data. Additional make sure you have enough RAM (>= 16 GB) so that running the pipeline doesn't take multiple days. A CPU that can run at least four threads will work with the code.

## The Folders

There are three folders in this project: the code used to generate figures; the code used for pipelines in both bash and Snakemake; and the tools used to assess pipeline recapitulation of ChIP-seq.

### Code to Generate Figures

This code holds unit-tested code for figure generation. Most likely the folder of least utility.

### Pipes

This folder contains all the pipeline code used in the paper, including the bash code used to test all 4560 pipelines as well as the Snakemake alternative for final pipelines. Interested users are provided with both, to use depending on their personal preference.

### ROCtool

This folder contains the code necessary to test a set of footprints against known ChIP-seq data. This ChIP-seq data is specific to GM12878; a new directory of ChIP-seq data needs to be made if a user produces footprint data from a different pipeline.