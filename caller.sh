#!/usr/bin/env bash

module purge
module load mash 
module load kraken2 
module load prokka 
module load nextflow 
module load signalp 
module load BLAST
module load  multiqc
module load panaroo

nextflow run main.nf -profile slurm #-resume
