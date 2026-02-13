rule sawfish_discover:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
    output:
        odir=directory(pjoin(WD, "{ref}", "sawfish-s{s}-q{q}")),
    log:
        time=pjoin(WD, "times", "{ref}", "sawfish1-s{s}-q{q}.time"),
    conda:
        "../envs/sawfish.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} sawfish discover --min-indel-size 50 --min-sv-mapq {wildcards.q} --threads {threads} --ref {input.fa} --bam {input.bam} --output-dir {output.odir}
        """


# --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XX.bed \
# --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz
rule sawfish_jointcall:
    input:
        odir=rules.sawfish_discover.output.odir,
    output:
        vcf=pjoin(WD, "{ref}", "sawfish-s{s}-q{q}.final", "genotyped.sv.vcf.gz"),
    params:
        odir=pjoin(WD, "{ref}", "sawfish-s{s}-q{q}.final"),
    log:
        time=pjoin(WD, "times", "{ref}", "sawfish2-s{s}-q{q}.time"),
    conda:
        "../envs/sawfish.yml"
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} sawfish joint-call --threads {threads} --min-sv-mapq {wildcards.q} --sample {input.odir} --output-dir {params.odir}
        """


rule sawfish_post:
    input:
        vcf=rules.sawfish_jointcall.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "sawfish-s{s}-q{q}.vcf.gz"),
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -Oz --include "INFO/SVTYPE='DEL' | INFO/SVTYPE='INS' | INFO/SVTYPE='DUP'" {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """
