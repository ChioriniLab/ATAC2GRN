# Pipelines README

This is the README for the pipelines used in this paper.

## Genomes Folder

Before using these pipelines, produce a genomes/hg19 or genomes/mm9 folder with a Bowtie2Index folder. Download [fasta files for the hg19 genome](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/) and concatenate them. Here is an option for command line input.

```
mkdir Bowtie2Index; cd Bowtie2Index
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
gzip -d *.fa.gz
cat *.fa > hg19.fasta
rm *.fa
bowtie2-build -f hg19.fasta hg19
rm hg19.fasta
```

Alternatively, download the genomes folder from the [NIDCR Box](https://nidcr.app.box.com/folder/0).

## Bash Pipeline

There is one final pipeline produced from this paper. That pipeline is the pipeline with the highest recovery of known ChIP-seq data (highest mean AUC).

With raw data files in data/your_data_name_here.1.fastq.gz and data/your_data_name_here.2.fastq.gz, and genome folders genomes/hg19 or genomes/mm9, the paired-end pipeline can be called with the command format

```
sh pipeline-v2-PE.sh -e name_of_output_folder_here -g hg19
```

Alternatively, with single-end read data in the format data/your_data_name_here.fastq.gz, the pipeline can be called as follows:

```
sh pipeline-v2-PE.sh -e name_of_output_folder_here -g hg19
```

If you have a CPU that supports fewer than 4 threads, some steps will be noticeably reduced in speed and you may want to adjust the code. Alternatively, if for some reason you can support many more than 4 threads, feel free to adjust the code to accomodate that.

## Snakemake Pipelines

If you are working on a cluster that uses slurm, the Snakefiles for these pipelines may be significantly more convenient. Snakemake allows pipelines to be broken up into many smaller steps, letting you use fewer cluster resources on simple tasks and making you a better neighbor on the cluster. Using Snakemake, you may need to adjust parameters and the cluster-v2.json file to the hardware of your cluster. Here at the NIH, we used the following format for the code:

```
module load snakemake
snakemake -s Snakefile-v2-hg19PE --unlock
snakemake -s Snakefile-v2-hg19PE -j56 -p -r --cluster "sbatch --time={cluster.time} --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus}" --cluster-config cluster-v2.json --keep-going --latency-wait 120 --rerun-incomplete
```

Both the VDI and the Singularity container have Snakemake, so you can load a Singularity image onto your cluster and use these pipelines.

# The Aim of These Pipelines

The point of these pipelines is to be usable out-of-the-box. Running the same code on two machines can be more difficult than it should be: if you encounter difficulties implementing these pipelines, please contact the corresponding author at jchiorini[at]dir[dot]nidcr[dot]nih[dot]gov or the first author at thomas[dot]pranzatelli[at]nih[dot]gov.