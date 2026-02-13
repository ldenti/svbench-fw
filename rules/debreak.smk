rule debreak:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
    output:
        vcf=pjoin(WD, "{ref}", "debreak-s{s}-q{q}", "debreak.vcf"),
    params:
        odir=pjoin(WD, "{ref}", "debreak-s{s}-q{q}"),
        s=lambda wildcards: f"--min_support {wildcards.s}" if wildcards.s != "0" else "",
    log:
        time=pjoin(WD, "times", "{ref}", "debreak-s{s}-q{q}.time"),
    conda:
        "../envs/debreak.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} debreak -t {threads} --ref {input.fa} --bam {input.bam} -o {params.odir} --min_size 50 --min_quality {wildcards.q} {params.s} --aligner minimap2
        """


rule debreak_post:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        vcf=rules.debreak.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "debreak-s{s}-q{q}.vcf.gz"),
    log:
        time=pjoin(WD, "times", "{ref}", "debreak-s{s}-q{q}.hiphase.time"),
    conda:
        "../envs/hiphase.yml"
    threads: workflow.cores
    shell:
        """
        echo {SAMPLE_NAME} > {input.vcf}.sample.txt
        bcftools reheader --samples {input.vcf}.sample.txt {input.vcf} | bcftools +fill-from-fasta -- -c REF -f {input.fa} | python3 ./scripts/to_upper.py | bgzip -c > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        /usr/bin/time -vo {log.time} hiphase --bam {input.bam} --reference {input.fa} --vcf {input.vcf}.gz --output-vcf {output.vcf} --threads {threads}
        """
