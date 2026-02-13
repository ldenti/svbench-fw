# rule svdss2_index:
#     input:
#         fa=REF,
#     output:
#         fmd=REF + ".fmd",
#     threads: workflow.cores
#     log:
#         time=pjoin(WD, "times", "rb3-index.time"),
#     conda:
#         "../envs/svdss2.yml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} SVDSS index -t{threads} -d {input.fa} -o {output.fmd}
#         """


rule svdss2_smooth:
    input:
        fa=REF,
        bam=BAM,  # pjoin(WD, "alignments.bam"),
    output:
        bam=pjoin(WD, "SVDSS2", "smoothed-q{q}.bam"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "svdss2-smooth-q{q}.time"),
    conda:
        "../envs/svdss2.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS smooth --reference {input.fa} --bam {input.bam} --threads {threads} --accp {wildcards.q} > {output.bam}
        samtools index {output.bam}
        """


rule svdss2_search:
    input:
        fmd=REF + ".fmd",
        bam=rules.svdss2_smooth.output.bam,
    output:
        sfs=pjoin(WD, "SVDSS2", "specifics-q{q}.tsv"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "svdss2-search-q{q}.time"),
    conda:
        "../envs/svdss2.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS search --index {input.fmd} --bam {input.bam} --threads {threads} > {output.sfs}
        """


rule svdss2_call:
    input:
        fa=REF,
        bam=pjoin(WD, "SVDSS2", "smoothed-q{q}.bam"),
        sfs=pjoin(WD, "SVDSS2", "specifics-q{q}.tsv"),
    output:
        vcf=pjoin(WD, "SVDSS2", "variations-q{q}-w{w}.vcf.gz"),
        sam=pjoin(WD, "SVDSS2", "variations-q{q}-w{w}.poa.sam"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "svdss2-call-q{q}-w{w}.time"),
    conda:
        "../envs/svdss2.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS call --min-sv-length 50 --min-cluster-weight {wildcards.w} --reference {input.fa} --bam {input.bam} --sfs {input.sfs} --threads {threads} --poa {output.sam} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule svdss2_gt:
    input:
        fa=REF,
        bam=BAM,
        vcf=rules.svdss2_call.output.vcf,
    output:
        vcf=pjoin(WD, "callsets", "SVDSS2-q{q}-w{w}.vcf.gz"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "SVDSS2-gt-q{q}-w{w}.time"),
    conda:
        "../envs/svdss2.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} kanpig gt --threads {threads} --input {input.vcf} --reads {input.bam} --reference {input.fa} | sed "s/FT/KF/g" | bcftools sort -Oz > {output.vcf}
        tabix -p vcf {output.vcf}
        """
