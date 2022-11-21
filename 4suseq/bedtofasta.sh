#!/bin/bash
#get sequences for bed file regions

indir=/scratch/Users/nila7826/4suseq/chip/fasta
outdir=/scratch/Users/nila7826/4suseq/chip/fasta
genomedir=/scratch/Users/nila7826/genome/hg38
genome=GRCh38.primary_assembly.genome.fa
bedfile=$1
name=$2
#software versions:
#bedtools v2.30.0
#-name option does not work on v2.28.0

bedtools getfasta \
-fo ${outdir}/${name}.fasta \
-name \
-fi ${genomedir}/${genome} \
-bed ${indir}/${bedfile}