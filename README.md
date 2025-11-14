# SVBench-FW

Modular and extensible framework to evaluate SV callers (from long reads) against SV truthsets created from diploid genomes.

The goal of this repository is to simplify SV calling performance evaluation while enabling a fair and standardised evaluation.

This framework is provided as a Snakemake pipeline. Additional python scripts can be used to summarize and plot the results.

### Prerequisites and usage
* conda/mamba
* snakemake
* seaborn (for optional plots)

``` sh
mamba create -c bioconda -c conda-forge -n svbench snakemake-minimal seaborn
conda activate svbench
# edit config/config.yml
snakemake -c 16 --use-conda --configfile config/config.yml -p [-n]

WD=$(grep "wd:" ./config/config.yml | cut -f2 -d" ")
ls $WD/*.csv
```

### Example
To test `svbench-fw`, we provide example data (zenodo) and the following instructions:
```
git clone https://github.com/ldenti/svbench-fw.git
cd svbench-fw

# mamba create -c bioconda -c conda-forge -n svbench snakemake-minimal
conda activate svbench

mkdir svbench-example
cd svbench-example
wget https://zenodo.org/records/17608271/files/svbench-fw.exampledata.tar.gz
tar xvfz svbench-fw.exampledata.tar.gz
bash write_config.sh > config.yml
cd ..
snakemake -c 4 --use-conda --configfile ./svbench-example/config.yml -p [-n]
ls ./svbench-example/SMK_OUT/*.csv
```

*Note 1:* this should take ~half an hour (using 4 threads)

*Note 2:* recall of all tools will be low since reads cover a small region of the chromosome whereas truthsets and contigs cover the entire chromosome.

### Supported tools
Read alignment:
* minimap2 (v2.28)

Small variants caller:
* deepvariant (v1.8.0)

Small variants phasing:
* whatshap (v2.5)

SV calling from long reads:
* SVDSS (v1.0.5, v2.1.0)
* sniffles (v2.3)
* cuteSV (v2.1.2)
* debreak (v1.3)
* SVision-pro (v2.4)
* Severus (v1.4.0)
* sawfish (v2.0.0)

SV calling from diploid assemblies:
* dipcall (v0.3)
* svim-asm (v1.0.3)
* hapdiff (commit e0abbb9a8095c70a0de23c49408a530901361b12)

Benchmarkers:
* truvari (v5.2.0)
* minda (commit 47d0fb5484b2b15865a94a9ba81436beaf52cf16)
  * no analysis on minda results

### Experiments
To replicate the experiments presented in the manuscript, follow the [instructions](https://github.com/ldenti/svbench/blob/main/plots/README.md) in the `plots` folder.
