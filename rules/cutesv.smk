rule cutesv:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        bed=pjoin(WD, "input", "trfs", "{ref}.bed"),
    output:
        vcf=pjoin(WD, "{ref}", "cutesv-s{s}-q{q}", "cutesv.vcf"),
    params:
        tmp=pjoin(WD, "{ref}", "cutesv-s{s}-q{q}.tmp"),
    log:
        time=pjoin(WD, "times", "{ref}", "cutesv-s{s}-q{q}.time"),
    conda:
        "../envs/cutesv.yml"
    threads: workflow.cores
    shell:
        """
        mkdir -p {params.tmp}
        /usr/bin/time -vo {log.time} cuteSV --sample {SAMPLE_NAME} --genotype --min_size 50 --threads {threads} --min_mapq {wildcards.q} --min_support {wildcards.s} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {input.fa} {output.vcf} {params.tmp}
        """


rule cutesv_post:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        vcf=rules.cutesv.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "cutesv-s{s}-q{q}.vcf.gz"),
    log:
        time=pjoin(WD, "times", "{ref}", "cutesv-s{s}-q{q}.hiphase.time"),
    conda:
        "../envs/hiphase.yml"
    threads: workflow.cores
    shell:
        """
        bcftools view --exclude "INFO/SVTYPE='BND'" {input.vcf} | python3 ./scripts/to_upper.py | bgzip -c > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        /usr/bin/time -vo {log.time} hiphase --bam {input.bam} --reference {input.fa} --vcf {input.vcf}.gz --output-vcf {output.vcf} --threads {threads}
        """
