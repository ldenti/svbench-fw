#!/bin/sh

SD=$(dirname $0)

set -xe

WD=$1

SAMPLE=NA12878 #HG002

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

    # for asm in giab hprc # this for HG002
    for asm in pp # this for NA12878
    do
	Rscript karyo.R $WD/$ref/asmcallsets-$asm/dipcall.bed $WD/$ref/asmcallsets-$asm/hapdiff.bed $ref $title $SAMPLE-$asm.$ref.pdf
    done
done
