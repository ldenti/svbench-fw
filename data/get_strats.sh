#!/bin/sh

WD=$1

pushd $WD

function gunzip_and_print {
    f=$1
    ref=$2

    bed=$(dirname $f)/$(basename $f .gz)

    if [ ! -f $bed ]
    then
        gunzip -kf $f
    fi
    bed=$(realpath $bed)

    echo "    $ref: $bed"
}

echo "vvvvvvv ADD THIS TO config.yml vvvvvvv"
echo "strat:"

echo "  hard:"
gunzip_and_print CHM13@all/Union/CHM13_alldifficultregions.bed.gz chm13
gunzip_and_print GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz hg38
gunzip_and_print GRCh37@all/Union/GRCh37_alldifficultregions.bed.gz hg19

echo "  easy:"
gunzip_and_print CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz chm13
gunzip_and_print GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz hg38
gunzip_and_print GRCh37@all/Union/GRCh37_notinalldifficultregions.bed.gz hg19

echo "  segdups:"
gunzip_and_print CHM13@all/SegmentalDuplications/CHM13_segdups.bed.gz chm13
gunzip_and_print GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz hg38
gunzip_and_print GRCh37@all/SegmentalDuplications/GRCh37_segdups.bed.gz hg19

echo "  lowmap:"
gunzip_and_print CHM13@all/Mappability/CHM13_lowmappabilityall.bed.gz chm13
gunzip_and_print GRCh38@all/Mappability/GRCh38_lowmappabilityall.bed.gz hg38
gunzip_and_print GRCh37@all/Mappability/GRCh37_lowmappabilityall.bed.gz hg19

echo "  lowcomplex:"
gunzip_and_print CHM13@all/LowComplexity/CHM13_AllTandemRepeatsandHomopolymers_slop5.bed.gz chm13
gunzip_and_print GRCh38@all/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz hg38
gunzip_and_print GRCh37@all/LowComplexity/GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed.gz hg19

echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

popd
