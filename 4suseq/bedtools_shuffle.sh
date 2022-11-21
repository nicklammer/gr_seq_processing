#!/bin/bash
#shuffle bed files (used to shuffle gene TSS bed file)

indir=/scratch/Users/nila7826/4suseq/chip/beds
outdir=/scratch/Users/nila7826/4suseq/chip/beds
genomedir=/scratch/Users/nila7826/genome/hg38
inbed=$1
name=$2
genomebed=hg38.genome
outname=${name}.tss.shuffle.bed
#software versions:
#bedtools v2.28.0

bedtools shuffle \
-i ${indir}/${inbed} \
-g ${genomedir}/${genomebed} |\
sort -T . -t $'\t' -k1,1 -k2,2n \
> ${outdir}/${outname}