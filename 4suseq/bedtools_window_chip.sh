#!/bin/bash
#bedtools window function to analyze overlap
#between gene TSSs and chip peaks
#return chip peaks with nearby TSS

genedir=/scratch/Users/nila7826/4suseq/chip/beds
chipdir=/scratch/Users/nila7826/4suseq/chip/beds
outdir=/scratch/Users/nila7826/4suseq/chip/fasta
genebed=$1
chipbed=$2
window=$3
outname=$4
#genebed is the gene TSS bed file
#chipbed is the chip peaks bed file
#window is the integer for window size
#software versions:
#bedtools v2.28.0

bedtools window \
-a ${chipdir}/${chipbed} \
-b ${genedir}/${genebed} \
-w ${window} \
-u \
> ${outdir}/${outname}