# Analyses

This folder contains all the scripts needed to reproduce the results presented in the manuscript. We assume you have run the workflow as described in the paper: individual HG002 with assemblies from giab and hprc and NA12878 with assemblies from platinum pedigree (pp).

All scripts/commands assume several tools and library to be installed (e.g., `truvari`, `seaborn`, `bedtools`, `R`...). Please refer also to `helper-asm.sh`, `helper-karyo.sh`, and `helper-longreads.sh` scripts.

All scripts can be run on the Snakemake working directory of both individuals (i.e., `$WD`), if not differently stated (e.g., `$HG002_WD`).

#### Assembly-based
``` sh
# Figure on assembly-based callsets
python3 ./plot_asm.py --asm giab [--confident] $WD

# Heatmap with pairwise similarity
python3 ./plot_comparison_asm.py $WD --asm giab [--confident] [--refine]

# Karyoplots
Rscript ./karyo.R $WD/chm13/asmcallsets-giab/dipcall.bed $WD/chm13/asmcallsets-giab/hapdiff.bed chm13 CHM13-giab karyo.pdf

# Chromosome with colored variants
python3 ./plot_chromosomes.py $HG002_WD/chm13/asmcallsets-giab/ chr3:143462850-143463100
```

#### Read-based vs assembly-based
``` sh
# stripplot
python3 ./plot_rank.py strip -a giab $HG002_WD/truvari.csv

# F1-based rankmap
python3 ./plot_rank.py heat -r chm13 -a giab $HG002_WD/truvari.csv

# delta(harmonized, notharmonized)
python3 ./plot_rank.py st $HG002_WD/truvari.csv

# F1 results on GIAB stratification
python3 ./plot_giabstrat.py -a giab --refine $HG002_WD/truvari.csv

# Parameters evalutions
python3 ./evaluate_parameters.py -a giab $HG002_WD/truvari.csv

# Singleton/multi stratification
python3 ./check_proximity.py csv $HG002_WD > proximity-refined.csv
python3 ./check_proximity.py pr proximity-refined.csv proximity-refined.pdf
```

#### Evaluation against GIAB (v5.0q and v0.6)
``` sh
# Statistics (#variants and length distribution)
python3 ./plot_giab_stats.py $HG002_WD

# P/R/F1 plot
python3 ./plot_giab.py [--refine] [--confident] $HG002_WD/truvari-giab.csv
```

#### Syntetic analysis
```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed

WD=HG002-SMK-WD # or NA12878-SMK-WD
for asm in giab hprc # pp for NA12878
do
  for caller in dipcall svim-asm hapdiff
  do
    bedtools intersect -header -a $WD/chm13/asmcallsets-$asm/$caller.vcf.gz -b chm13v2-unique_to_hg38.bed -u > $WD/chm13/asmcallsets-$asm/$caller.non-syntenic.vcf
    bedtools intersect -header -a $WD/chm13/asmcallsets-$asm.confident/$caller.vcf.gz -b chm13v2-unique_to_hg38.bed -u > $WD/chm13/asmcallsets-$asm.confident/$caller.non-syntenic.vcf
  done
done

python3 ./plot_asm_syntenic.py --asm giab [--confident] $HG002-SMK-WD # or NA12878-SMK-WD 

# Callers F1 comparison on syntetic regions, CHM13, GRCh38
bash ./evaluate_nonsyntetic.sh $WD chm13v2-unique_to_hg38.bed
```

#### Assembly-based SVs per GIAB tier
```
# Table summarizing stratifications
python3 analyze_strats.py [HG002-SMK-WD]

# Plot with SV length per stratifications (CHM13, giab)
WD=HG002-SMK-WD # or NA12878-SMK-WD
for bed in $WD/input/strats/*/chm13.bed
do
  region=$(basename $(dirname $bed))
  for asm in giab hprc # pp for NA12878
  do
    for caller in dipcall svim-asm hapdiff
    do
      bedtools intersect -header -a $WD/chm13/asmcallsets-$asm/$caller.vcf.gz -b $bed -u > $WD/chm13/asmcallsets-$asm/$caller.$region.vcf
      bedtools intersect -header -a $WD/chm13/asmcallsets-$asm.confident/$caller.vcf.gz -b $bed -u > $WD/chm13/asmcallsets-$asm.confident/$caller.$region.vcf
    done
  done
done

python3 plot_asm_strats.py --asm giab [--confident] $WD
```
