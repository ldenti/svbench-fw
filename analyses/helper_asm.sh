#!/bin/sh

SD=$(dirname $0)

set -xe

HG=$1
NA=$2

# # ASM-BASED BIG FIGURE
# #######################

# for asm in giab hprc
# do
#     python3 $SD/plot_asm.py $HG --asm $asm             -o hg002-$asm.full.pdf
#     python3 $SD/plot_asm.py $HG --asm $asm --confident -o hg002-$asm.conf.pdf
# done

# asm=pp
# python3 $SD/analyses/plot_asm.py $NA --asm $asm             -o na12878-$asm.full.pdf
# python3 $SD/plot_asm.py $NA --asm $asm --confident -o na12878-$asm.conf.pdf


# # ASM-BASED STRATIFICATION
# ###########################

# for asm in giab hprc
# do
#     python3 $SD/plot_asm_strats.py $HG --asm $asm             -o strats.hg002-$asm-full.pdf
#     python3 $SD/plot_asm_strats.py $HG --asm $asm --confident -o strats.hg002-$asm-conf.pdf
# done

# asm=pp
# python3 $SD/plot_asm_strats.py $NA --asm $asm             -o strats.na12878-$asm-full.pdf
# python3 $SD/plot_asm_strats.py $NA --asm $asm --confident -o strats.na12878-$asm-conf.pdf

# ######################################################################

# # ASM-BASED NON-SYNTENIC
# #########################

# for asm in giab hprc
# do
#     python3 $SD/plot_asm_syntenic.py $HG --asm $asm             -o syntenic.hg002-$asm-full.pdf
#     python3 $SD/plot_asm_syntenic.py $HG --asm $asm --confident -o syntenic.hg002-$asm-conf.pdf
# done

# asm=pp
# python3 $SD/plot_asm_syntenic.py $NA --asm $asm             -o syntenic.na12878-$asm-full.pdf
# python3 $SD/plot_asm_syntenic.py $NA --asm $asm --confident -o syntenic.na12878-$asm-conf.pdf

# ######################################################################

# ASM-BASED COMPARISON HEATMAPS
################################

for asm in giab hprc
do
    python3 $SD/plot_comparison_asm.py $HG                      --asm $asm -o heatmap.hg002-${asm}.full-noref.pdf
    python3 $SD/plot_comparison_asm.py $HG             --refine --asm $asm -o heatmap.hg002-${asm}.full-ref.pdf
    python3 $SD/plot_comparison_asm.py $HG --confident          --asm $asm -o heatmap.hg002-${asm}.conf-noref.pdf
    python3 $SD/plot_comparison_asm.py $HG --confident --refine --asm $asm -o heatmap.hg002-${asm}.conf-ref.pdf
done

asm=pp
python3 $SD/plot_comparison_asm.py $NA                      --asm $asm -o heatmap.na12878-${asm}.full-noref.pdf
python3 $SD/plot_comparison_asm.py $NA             --refine --asm $asm -o heatmap.na12878-${asm}.full-ref.pdf
python3 $SD/plot_comparison_asm.py $NA --confident          --asm $asm -o heatmap.na12878-${asm}.conf-noref.pdf
python3 $SD/plot_comparison_asm.py $NA --confident --refine --asm $asm -o heatmap.na12878-${asm}.conf-ref.pdf

######################################################################

# # ASM-BASED CONFIDENT REGIONS
# ##############################

# for ref in chm13 hg38 hg19
# do
#     for asm in giab hprc
#     do
# 	Rscript $SD/karyo.R $HG/$ref/asmcallsets-$asm/dipcall.bed $HG/$ref/asmcallsets-$asm/hapdiff.bed $ref $ref-$asm HG002-$asm.$ref.pdf
#     done
# done

# for ref in chm13 hg38 hg19
# do
#     for asm in pp
#     do
# 	Rscript $SD/karyo.R $NA/$ref/asmcallsets-$asm/dipcall.bed $NA/$ref/asmcallsets-$asm/hapdiff.bed $ref $ref-$asm NA12878-$asm.$ref.pdf
#     done
# done
