#!/bin/sh

SD=$(dirname $0)

set -xe

WD=$1 # SNAKEMAKE DIRECTORY

SAMPLE=NA12878 # HG002
ASM=(pp) # (giab hprc)

for ref in chm13 hg38 hg19
do
    title=""
    if [[ $ref == "chm13" ]]
    then
	title="T2T-CHM13"
    elif [[ $ref == "hg38" ]]
    then
	title="GRCh38"
    elif [[ $ref == "hg19" ]]
    then
	title="GRCh37"
    fi

    for asm in "${ASM[@]}"
    do
	Rscript $SD/karyo.R $WD/$ref/asmcallsets-$asm/dipcall.bed $WD/$ref/asmcallsets-$asm/hapdiff.bed $ref $title $SAMPLE-$asm.$ref.pdf
    done
done
