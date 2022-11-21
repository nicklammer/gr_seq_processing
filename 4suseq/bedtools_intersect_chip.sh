#!/bin/bash

indir=/scratch/Users/nila7826/chip_stallcup
outdir=/scratch/Users/nila7826/chip_stallcup
bed1=$1
bed2=$2
name=$3

#software versions:
#bedtools v2.28.0

bedtools intersect \
-a ${indir}/${bed1} \
-b ${indir}/${bed2} \
-sorted \
> ${outdir}/${name}.intersect.bed