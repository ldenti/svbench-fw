
rule severus:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        bed=pjoin(WD, "input", "trfs", "{ref}.bed"),
        vcf=pjoin(WD, "{ref}", "deepvariant-phased.vcf.gz"),
    output:
        vcf=pjoin(WD, "{ref}", "severus-s{s}-q{q}", "all_SVs", "severus_all.vcf"),
    params:
        tmp=pjoin(WD, "{ref}", "severus-s{s}-q{q}"),
    log:
        time=pjoin(WD, "times", "{ref}", "severus-s{s}-q{q}.time"),
    conda:
        "../envs/severus.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} severus --target-sample {SAMPLE_NAME} --min-sv-size 50 --min-mapq {wildcards.q} --min-support {wildcards.s} --target-bam {input.bam} --vntr-bed {input.bed} --out-dir {params.tmp} -t {threads} --phasing-vcf {input.vcf}
        """


rule severus_post:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        vcf=rules.severus.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "severus-s{s}-q{q}.vcf.gz"),
    log:
        time=pjoin(WD, "times", "{ref}", "severus-s{s}-q{q}.hiphase.time"),
    conda:
        "../envs/hiphase.yml"
    threads: workflow.cores
    shell:
        """
        bcftools view --exclude "INFO/SVTYPE='BND'" {input.vcf} | bcftools +fill-from-fasta -- -c REF -f {input.fa} | python3 ./scripts/to_upper.py | sed "s/minimap2/{SAMPLE_NAME}/g" | bgzip -c > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        /usr/bin/time -vo {log.time} hiphase --bam {input.bam} --reference {input.fa} --vcf {input.vcf}.gz --output-vcf {output.vcf} --threads {threads}
        """
