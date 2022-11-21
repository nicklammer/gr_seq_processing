#!/bin/bash

workdir=/scratch/Users/nila7826/chip_stallcup/nfcore
config=$1
paramfile=$2
#config contains SLURM-specific parameters
#paramfile is a json file with nf-core parameters
#software versions:
#singularity v3.1.1

/Users/nila7826/nextflow/nextflow run nf-core/chipseq -r 1.2.2 \
-profile singularity \
-c ${config} \
-work-dir ${workdir} \
-params-file ${paramfile}