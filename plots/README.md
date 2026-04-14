# Analyses

This folder contains all the scripts needed to reproduce the results presented in the manuscript. All the scripts assume `truvari` and `seaborn` to be installed (e.g., via conda).

We ran the Snakemake pipeline three times, one per reference genome (i.e., hg19, hg38, t2t). For each reference genome, update the [configuration `config/config.yml`](https://github.com/ldenti/svbench/blob/main/config/config.yml) and run the Snakemake pipeline. We assume the variables `{hg19,hg38,t2t}_smk_wd` to be the path to the output directories of the three Snakemake runs (one per reference sequence).

#### Assembly-based analysis
The folder `truth` in the snakemake output directory contains the assembly-based truth sets. They can be analyzed using the following commands:
``` sh
# Figure on assembly-based callsets
python3 ./plot_asm.py [--confident] $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd

# Three genomes heatmap (jaccard similarity/accuracy)
python3 ./plot_comparison_asm.py [--refine] $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
```

#### Read-based vs assembly-based
``` sh
# summarize F1 in a single plot
python3 ./plot_truvari.py all $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
# rank map depending on F1
python3 ./plot_truvari.py rank $t2t_smk_wd T2T

# F1 results on GIAB stratification
python3 ./plot_giabstrat.py [--refine] $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
```

#### Evaluation against GIAB (v1.1 and v0.6)
``` sh
python3 ./plot_giab.py [--refine] [--bed] $t2t_smk_wd $hg38_smk_wd $hg19_smk_wd
```

<!--
## Double assembly analyses
These scripts analyze the results obtained from both assemblies:
1. run the snakemake on 3 references using HPRC contigs
2. run the snakemake on 3 references using GIAB haplotypes (you can symlink the callsets directory)
3. run the GIAB v1.1 scripts (see above) from one of the two runs (since this is independent from the assembly used)

#### Pairwise comparison of truth sets
This will compare all assembly-based truth sets and GIAB v1.1 and v0.6 (on hg19).
``` sh
bash ./run_truvari.sh t2t.output_directory t2t.smk_workdir_on_hprc t2t.smk_workdir_on_giab hg38.giab-v11.vcf.gz .
bash ./run_truvari.sh hg38.output_directory hg38.smk_workdir_on_hprc hg38.smk_workdir_on_giab hg38.giab-v11.vcf.gz .
bash ./run_truvari.sh hg19.output_directory hg19.smk_workdir_on_hprc hg19.smk_workdir_on_giab hg19.giab-v11.vcf.gz hg19.giab-v06.vcf.gz

# plot the accuracy heatmap
python3 ./plot_comparison_asm_full.py t2t.output_directory hg38.output_directory hg19.output_directory
```

#### Full rankmap
This will produce a rankmap containing both assemblies and the GIAB v1.1 truth set.
``` sh
python3 scripts/plot_rankmap.py hg38.smk_workdir_on_hprc hg38.smk_workdir_on_giab PlotTitle
# ^ adapt for other references
```
-->

