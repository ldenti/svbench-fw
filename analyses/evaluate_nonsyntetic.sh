#!/bin/sh

set -e

WD=$1
synt_bed=$2

OD=$WD/chm13/synthetic_evaluation

ref=$WD/input/refs/chm13.fa

sort -k1,1 -k2,2n $synt_bed > $synt_bed.sorted
sort -k1,1 $ref.fai > $WD/chm13.sorted.fai

bedtools complement -i $synt_bed.sorted -g $WD/chm13.sorted.fai > $synt_bed.complement.bed

for truth in dipcall svim-asm hapdiff
do
    mkdir -p $OD/$truth

    truth_vcf=$WD/chm13/asmcallsets-giab/$truth.vcf.gz
    for caller in cutesv-s4-q20 debreak-s0-q20 sawfish-s0-q20 severus-s4-q20 sniffles-s0-q20 SVDSS-s4-q0 svisionpro-s4-q20
    do
	if [ ! -d $OD/$truth/$caller ]
	then
	    caller_vcf=$WD/chm13/callsets/$caller.vcf.gz
	    truvari bench --includebed $synt_bed.complement.bed --passonly --pick ac --dup-to-ins --reference $ref --base $truth_vcf --comp $caller_vcf --output $OD/$truth/$caller
            truvari refine --reference $ref --coords R --use-original-vcfs --threads 16 --align mafft $OD/$truth/$caller
            truvari ga4gh --input $OD/$truth/$caller --output $OD/$truth/$caller/ga4gh_with_refine
	fi
    done
done

echo "Truth,Caller,P_full,R_full,F1_full,P_synt,R_synt,F1_synt,P_GRCh38,R_GRCh38,F1_GRCh38"
for truth in dipcall svim-asm hapdiff
do
    for caller in cutesv-s4-q20 debreak-s0-q20 sawfish-s0-q20 severus-s4-q20 sniffles-s0-q20 SVDSS-s4-q0 svisionpro-s4-q20
    do
	fp=$WD/chm13/truvari/giab/$truth/$caller/ga4gh_with_refine.summary.json
	p1=$(grep "precision" $fp | cut -f2 -d":" | cut -c2-6)
	r1=$(grep "recall" $fp | cut -f2 -d":" | cut -c2-6)
	f1=$(grep "f1" $fp | cut -f2 -d":" | cut -c2-6)

	fp=$OD/$truth/$caller/ga4gh_with_refine.summary.json
	p2=$(grep "precision" $fp | cut -f2 -d":" | cut -c2-6)
	r2=$(grep "recall" $fp | cut -f2 -d":" | cut -c2-6)
	f2=$(grep "f1" $fp | cut -f2 -d":" | cut -c2-6)

	# dp=$(printf "%0.2f\n" "$(echo "$p2-$p1" | bc -l)")
	# dr=$(printf "%0.2f\n" "$(echo "$r2-$r1" | bc -l)")
	# df=$(printf "%0.2f\n" "$(echo "$f2-$f1" | bc -l)")

	fp=$WD/hg38/truvari/giab/$truth/$caller/ga4gh_with_refine.summary.json
	p3=$(grep "precision" $fp | cut -f2 -d":" | cut -c2-6)
	r3=$(grep "recall" $fp | cut -f2 -d":" | cut -c2-6)
	f3=$(grep "f1" $fp | cut -f2 -d":" | cut -c2-6)

	echo $truth,$caller,$p1,$r1,$f1,$p2,$r2,$f2,$p3,$r3,$f3
    done
done
