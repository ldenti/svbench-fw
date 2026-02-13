rule svdss2_index:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
    output:
        fmd=pjoin(WD, "input", "refs", "{ref}.fa.fmd"),
    threads: workflow.cores
    log:
        time=pjoin(WD, "times", "{ref}", "rb3-index.time"),
    conda:
        "../envs/svdss2.yml"
    shell:
        """
        /usr/bin/time -vo {log.time} SVDSS index -t{threads} -d {input.fa} -o {output.fmd}
        """


# rule svdss2_smooth_ht:
#     input:
#         fa=REF,
#         bam=BAM_HT,
#     output:
#         bam=pjoin(WD, "SVDSS2ht", "smoothed-q{q}.bam"),
#     threads: workflow.cores
#     log:
#         time=pjoin(WD, "times", "svdss2ht-smooth-q{q}.time"),
#     conda:
#         "../envs/svdss2.yml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} SVDSS smooth --reference {input.fa} --bam {input.bam} --threads {threads} --accp {wildcards.q} > {output.bam}
#         samtools index {output.bam}
#         """
# rule svdss2_search_ht:
#     input:
#         fmd=REF + ".fmd",
#         bam=rules.svdss2_smooth_ht.output.bam,
#     output:
#         sfs=pjoin(WD, "SVDSS2ht", "specifics-q{q}.tsv"),
#     threads: workflow.cores
#     log:
#         time=pjoin(WD, "times", "svdss2ht-search-q{q}.time"),
#     conda:
#         "../envs/svdss2.yml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} SVDSS search --index {input.fmd} --bam {input.bam} --threads {threads} > {output.sfs}
#         """
# rule svdss2_call_ht:
#     input:
#         fa=REF,
#         bam=pjoin(WD, "SVDSS2ht", "smoothed-q{q}.bam"),
#         sfs=pjoin(WD, "SVDSS2ht", "specifics-q{q}.tsv"),
#     output:
#         vcf=pjoin(WD, "SVDSS2ht", "variations-q{q}-w{w}.vcf.gz"),
#         sam=pjoin(WD, "SVDSS2ht", "variations-q{q}-w{w}.poa.sam"),
#     threads: workflow.cores
#     log:
#         time=pjoin(WD, "times", "svdss2-call-q{q}-w{w}.time"),
#     conda:
#         "../envs/svdss2.yml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} SVDSS call --min-sv-length 50 --min-cluster-weight {wildcards.w} --reference {input.fa} --bam {input.bam} --sfs {input.sfs} --threads {threads} --poa {output.sam} | bgzip -c > {output.vcf}
#         tabix -p vcf {output.vcf}
#         """
# rule svdss2_gt_ht:
#     input:
#         fa=REF,
#         bam=BAM_HT,
#         vcf=rules.svdss2_call_ht.output.vcf,
#     output:
#         vcf=pjoin(WD, "callsets", "SVDSS2ht-q{q}-w{w}.vcf.gz"),
#     threads: workflow.cores
#     log:
#         time=pjoin(WD, "times", "SVDSS2ht-gt-q{q}-w{w}.time"),
#     conda:
#         "../envs/svdss2.yml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} kanpig gt --threads {threads} --input {input.vcf} --reads {input.bam} --reference {input.fa} | sed "s/FT/KF/g" | bcftools sort -Oz > {output.vcf}
#         tabix -p vcf {output.vcf}
#         """
