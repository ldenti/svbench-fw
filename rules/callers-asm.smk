min_l = 50

# rule intersect_bed:
#     input:
#         bed1=DIPBED,
#         bed2=HAPBED,
#     output:
#         bed=BED,
#     conda:
#         "../envs/bedtools.yml"
#     shell:
#         """
#         bedtools intersect -a {input.bed1} -b {input.bed2} > {output.bed}
#         """


rule dipcall:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        # bed=PARBED,
        hap1=pjoin(WD, "input", "asm-{a}", "hap1.fa"),
        hap2=pjoin(WD, "input", "asm-{a}", "hap2.fa"),
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall", "prefix.dip.vcf.gz"),
        bed=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall", "prefix.dip.bed"),
        bam1=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall", "prefix.hap1.bam"),
        bam2=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall", "prefix.hap2.bam"),
    params:
        wdir=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall"),
        prefix=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall", "prefix"),
        mak=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.mak"),
    threads: workflow.cores
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        run-dipcall -t {threads} {params.prefix} {input.fa} {input.hap1} {input.hap2} > {params.mak}
        mkdir -p {params.wdir}
        make -j 2 -f {params.mak}
        tabix -p vcf {output.vcf}
        samtools index {output.bam1}
        samtools index {output.bam2}
        """


rule clean_dipcall:
    input:
        fai=pjoin(WD, "input", "refs", "{ref}.fa.fai"),
        bed=rules.dipcall.output.bed,
        vcf=rules.dipcall.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.vcf.gz"),
        bed=pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.bed"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools reheader --fai {input.fai} {input.vcf} | bcftools norm --multiallelics - | bcftools view -v indels -i '(ILEN <= -{min_l} || ILEN >= {min_l})' -Oz > {output.vcf}
        tabix -p vcf {output.vcf}
        cp {input.bed} {output.bed}
        """


######################################################################
######################################################################
######################################################################


rule svim_asm:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam1=rules.dipcall.output.bam1,
        bam2=rules.dipcall.output.bam2,
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "svim-asm", "variants.vcf"),
    params:
        wd=pjoin(WD, "{ref}", "asmcallsets-{a}", "svim-asm"),
    conda:
        "../envs/svimasm.yml"
    shell:
        """
        mkdir -p {params.wd}
        svim-asm diploid --min_sv_size {min_l} {params.wd} {input.bam1} {input.bam2} {input.fa}
        """


rule clean_svim_asm:
    input:
        vcf=rules.svim_asm.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "svim-asm.vcf.gz"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools view -Oz -v indels -i '(ILEN <= -{min_l} || ILEN >= {min_l})' {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule get_hapdiff:
    output:
        exe=pjoin(WD, "software", "hapdiff", "hapdiff.py"),
        repo=directory(pjoin(WD, "software", "hapdiff")),
    shell:
        """
        git clone https://github.com/KolmogorovLab/hapdiff {output.repo}
        cd {output.repo}
        git checkout e0abbb9a8095c70a0de23c49408a530901361b12
        git submodule update --init
        make
        """


rule hapdiff:
    input:
        exe=rules.get_hapdiff.output.exe,
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        hap1=pjoin(WD, "input", "asm-{a}", "hap1.fa"),
        hap2=pjoin(WD, "input", "asm-{a}", "hap2.fa"),
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff", "hapdiff_phased.vcf.gz"),
        bed=pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff", "confident_regions.bed"),
        bed2=pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff.bed"),
    params:
        outd=pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff"),
    conda:
        "../envs/hapdiff.yml"
    threads: workflow.cores
    shell:
        """
        {input.exe} --reference {input.fa} --pat {input.hap1} --mat {input.hap2} --out-dir {params.outd} -t {threads}
        cp {input.bed} {input.bed2}
        """


rule hapdiff_post:
    input:
        vcf=rules.hapdiff.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff.vcf.gz"),
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools view -Oz -v indels -i '(ILEN <= -{min_l} || ILEN >= {min_l})' {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


######################################################################
######################################################################
######################################################################


rule subset:
    input:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.vcf.gz"),
        bed=lambda wildcards: (
            rules.hapdiff.output.bed
            if wildcards.asmcaller == "hapdiff"
            else rules.dipcall.output.bed
        ),
    output:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}.confident", "{asmcaller}.vcf.gz"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -header -u -a {input.vcf} -b {input.bed} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
