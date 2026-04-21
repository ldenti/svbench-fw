#!/bin/sh

SD=$(dirname $0)

set -xe

HG=$1
NA=$2

python3 $SD/plot_rank.py strip -a giab $HG/truvari.csv -o hg002-giab.strip.pdf
python3 $SD/plot_rank.py strip -a hprc $HG/truvari.csv -o hg002-hprc.strip.pdf
python3 $SD/plot_rank.py strip -a pp   $NA/truvari.csv -o na12878-pp.strip.pdf

for ref in chm13 hg38 hg19
do
    python3 $SD/plot_rank.py heat -r $ref -a giab $HG/truvari.csv -o hg002-giab.rankmap-$ref.pdf
    python3 $SD/plot_rank.py heat -r $ref -a hprc $HG/truvari.csv -o hg002-hprc.rankmap-$ref.pdf
    python3 $SD/plot_rank.py heat -r $ref -a pp   $NA/truvari.csv -o na12878-pp.rankmap-$ref.pdf
done

python3 $SD/plot_giabstrat.py -a giab --refine $HG/truvari.csv -o hg002-giab.strat.pdf
python3 $SD/plot_giabstrat.py -a hprc --refine $HG/truvari.csv -o hg002-hprc.strat.pdf
# python3 $SD/plot_giabstrat.py -a pp   --refine $NA/truvari.csv

python3 $SD/plot_giab.py --refine             $HG/truvari-giab.csv -o full.pdf
python3 $SD/plot_giab.py --refine --confident $HG/truvari-giab.csv -o conf.pdf
