#!/bin/bash
#deduplicate aligned reads with UMI tools

indir=/scratch/Users/nila7826/4suseq/align
outdir=/scratch/Users/nila7826/4suseq/dedup
bam=$1
name=$2
#bam is the aligned read file as .bam
#software versions:
#python v3.7.4
#UMI-tools v1.1.2

umi_tools dedup \
-I ${indir}/${bam} \
--paired \
--unpaired-reads discard \
--chimeric-pairs discard \
-S ${outdir}/${name}.dedup.bam