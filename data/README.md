# Data

### Reference genomes and annotations
```
# GRCh37
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
samtools faidx human_g1k_v37.fasta
wget https://github.com/fritzsedlazeck/Sniffles/blob/master/annotations/human_hs37d5.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/genome-stratifications-GRCh37@all.tar.gz
tar xvfz genome-stratifications-GRCh37@all.tar.gz

# GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/genome-stratifications-GRCh38@all.tar.gz
tar xvfz genome-stratifications-GRCh38@all.tar.gz

# CHM13-T2T
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz
samtools faidx chm13v2.0.fa
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/genome-stratifications-CHM13@all.tar.gz
tar xvfz genome-stratifications-CHM13@all.tar.gz

# To unzip the stratifications we need and produce the corresponding config file
bash get_strats.sh .
```

### Assemblies
```
# HG002 (GIAB)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
# split fasta in maternal and paternal using samtools faidx

# HG002 (HPRC)
curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc
agc getset ../../HPRC-yr1.agc HG002.1 > HG002.1.fa
samtools faidx HG002.1.fa
agc getset ../../HPRC-yr1.agc HG002.2 > HG002.2.fa
samtools faidx HG002.2.fa

# NA12878 (platinum pedigree)
wget https://42basepairs.com/download/s3/platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype1.fasta{.fai}
wget https://42basepairs.com/download/s3/platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype2.fasta{.fai}
```

### HiFi samples
```
# HG002 (HPRC)
REPO=https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/
wget $REPO/m64012_190920_173625.Q20.fastq
wget $REPO/m64012_190921_234837.Q20.fastq
wget $REPO/m64015_190920_185703.Q20.fastq
wget $REPO/m64015_190922_010918.Q20.fastq
cat *.fastq | gzip -c > hprc-hifi2.0-15k.fq.gz

# NA12878 (platinum pedigree)
wget https://42basepairs.com/download/s3/platinum-pedigree-data/data/hifi/mapped/CHM13/NA12878.CHM13.haplotagged.bam{.bai}
samtools fastq NA12878.CHM13.haplotagged.bam | gzip -c > NA12878.CHM13.haplotagged.fq.gz
```

### GIAB curated SV callsets
```
# Get GIAB v0.6
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.{vcf.gz,vcf.gz.tbi,bed}

# Get GIAB v5.0
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_{GRCh37,GRCh38,CHM13v2.0}_v5.0q_stvar.vcf.gz{.tbi}
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_{GRCh37,GRCh38,CHM13v2.0}_v5.0q_stvar.benchmark.bed

# get SVs
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' HG002_CHM13v2.0_v5.0q_stvar.vcf.gz > chm13-giab5.sv.vcf.gz
tabix -p vcf chm13-giab5.sv.vcf.gz
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' HG002_GRCh38_v5.0q_stvar.vcf.gz > hg38-giab5.sv.vcf.gz
tabix -p vcf hg38-giab5.sv.vcf.gz
bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' HG002_GRCh37_v5.0q_stvar.vcf.gz > hg19-giab5.sv.vcf.gz
tabix -p vcf hg19-giab5.sv.vcf.gz
```

<!-- Original VCFs are provided by Severus people ([variant_calls_and_benchmarks.tar.gz](https://zenodo.org/records/14541057)). To run truvari bench, we had to remove BND from some callsets (sniffles2, cuteSV, and severus):
```
bcftools view -Oz --exclude "INFO/SVTYPE='BND'" input.vcf.gz > output.vcf.gz
tabix -p vcf output.vcf.gz
``` -->
