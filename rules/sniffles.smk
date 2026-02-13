rule sniffles:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        bed=pjoin(WD, "input", "trfs", "{ref}.bed"),
    output:
        vcf=pjoin(WD, "{ref}", "sniffles-s{s}-q{q}", "sniffles.vcf"),
    params:
        s=lambda wildcards: wildcards.s if wildcards.s != "0" else "auto",
    log:
        time=pjoin(WD, "times", "{ref}", "sniffles-s{s}-q{q}.time"),
    conda:
        "../envs/sniffles.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} sniffles --threads {threads} --mapq {wildcards.q} --minsupport {params.s} --phase --minsvlen 50 --reference {input.fa} --input {input.bam} --vcf {output.vcf} --tandem-repeats {input.bed}
        """


rule sniffles_post:
    input:
        vcf=rules.sniffles.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "sniffles-s{s}-q{q}.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --exclude "INFO/SVTYPE='BND'" {input.vcf} | python3 ./scripts/remove_strange_symbolic.py | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
