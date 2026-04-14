rule svdss_smooth:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
    output:
        bam=pjoin(WD, "{ref}", "SVDSS", "smoothed.selective.bam"),
    params:
        wd=pjoin(WD, "{ref}", "SVDSS"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "{ref}", "svdss-smooth.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS smooth --reference {input.fa} --bam {input.bam} --workdir {params.wd} --threads {threads}
        samtools index {output.bam}
        """


rule svdss_search:
    input:
        fmd=pjoin(WD, "input", "refs", "{ref}.fa.fmd"),
        bam=rules.svdss_smooth.output.bam,
    output:
        sfs=pjoin(WD, "{ref}", "SVDSS", "solution_batch_0.assembled.sfs"),
    params:
        wd=pjoin(WD, "{ref}", "SVDSS"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "{ref}", "svdss-search.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS search --index {input.fmd} --bam {input.bam} --threads {threads} --workdir {params.wd} --assemble
        """


rule svdss_call:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=rules.svdss_smooth.output.bam,
        sfs=rules.svdss_search.output.sfs,
    output:
        vcf=pjoin(WD, "{ref}", "SVDSS", "s{s}", "svs_poa.vcf"),
    params:
        wd=pjoin(WD, "{ref}", "SVDSS"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "{ref}", "svdss-call-s{s}-q0.time"),
    conda:
        "../envs/svdss.yml"
    shell:
        """
        n=$(ls {params.wd}/solution_batch_*.assembled.sfs | wc -l)
        # NOTE: filter on length is "wrong" (>) in SVDSS, need to set to 49
        /usr/bin/time -vo {log.time} SVDSS call --reference {input.fa} --bam {input.bam} --threads {threads} --workdir {params.wd} --batches ${{n}} --min-cluster-weight {wildcards.s} --min-sv-length 49
        cp {params.wd}/svs_poa.vcf {output.vcf}
        """


rule svdss_post:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=rules.svdss_smooth.output.bam,
        vcf=rules.svdss_call.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "SVDSS-s{s}-q0.vcf.gz"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "{ref}", "svdss-hiphase-s{s}-q0.time"),
    conda:
        "../envs/hiphase.yml"
    shell:
        """
        echo {SAMPLE_NAME} > {input.vcf}.sample.txt
        bcftools reheader --samples {input.vcf}.sample.txt {input.vcf} | python3 ./scripts/to_upper.py | bgzip -c > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        /usr/bin/time -vo {log.time} hiphase --bam {input.bam} --reference {input.fa} --vcf {input.vcf}.gz --output-vcf {output.vcf} --threads {threads} > {output.vcf}
        """
