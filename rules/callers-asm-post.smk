wildcard_constraints:
    h="(1)|(2)",
    mode="(full)|(conf)",
    confident="()|(.confident)",


# === Pairwise comparison
# ===========================
rule assemblybased_pairwise_truvari:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf1=pjoin(WD, "{ref}", "asmcallsets-{a}", "{truth1}.vcf.gz"),
        vcf2=pjoin(WD, "{ref}", "asmcallsets-{a}", "{truth2}.vcf.gz"),
        bed=lambda wildcards: (
            pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff", "confident_regions.bed")
            if wildcards.truth1 == "hapdiff"
            else pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.bed")
        ),
    output:
        d=directory(
            pjoin(
                WD,
                "{ref}",
                "asmcallsets-{a}",
                "comparison-{mode}",
                "{truth2}-against-{truth1}",
            )
        ),
        dd=directory(
            pjoin(
                WD,
                "{ref}",
                "asmcallsets-{a}",
                "comparison-{mode}",
                "{truth2}-against-{truth1}",
                "phab_bench",
            )
        ),
    params:
        bed=lambda wildcards, input: (
            "--includebed " + input.bed if wildcards.mode == "conf" else ""
        ),
    conda:
        "../envs/truvari.yml"
    threads: workflow.cores / 2
    shell:
        """
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins -s 50 -S 50 {params.bed} --reference {input.fa} --base {input.vcf1} --comp {input.vcf2} --output {output.d}
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        if [ -d {output.dd} ]
        then
            truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine
        else
            mkdir -p {output.dd}
            touch {output.d}/NOREGIONSTOREFINE
            cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
        fi
        """


# === TTmars-like analysis
# ============================
rule remove_info:
    input:
        vcf=pjoin(WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.vcf.gz"),
    output:
        vcf=pjoin(
            WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.noinfo.vcf.gz"
        ),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        bash ./scripts/remove_info_from_vcf.sh {input.vcf} | bcftools view -Oz -v indels -i '(ILEN <= -30 || ILEN >= 30)' > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vcf2regions:
    input:
        vcf=rules.remove_info.output.vcf,
    output:
        txt=pjoin(
            WD,
            "{ref}",
            "asmcallsets-{a}{confident}",
            "{asmcaller}.noinfo.regions-w{w}.txt",
        ),
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        python3 ./scripts/vcf2regions.py {input.vcf} {wildcards.w} > {output.txt}
        """


rule get_contigs:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=rules.remove_info.output.vcf,
        txt=rules.vcf2regions.output.txt,
    output:
        fa=pjoin(
            WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.hap{h}-w{w}.fa"
        ),
    params:
        tag=lambda wildcards: "first" if wildcards.h == "1" else "second",
    conda:
        "../envs/sambcftools.yml"
    shell:
        """
        grep {params.tag} {input.txt} | cut -f1 -d" " | while read region ; do samtools faidx {input.fa} $region | bcftools consensus {input.vcf} -H {wildcards.h} ; done > {output.fa} 2> {output.fa}.log
        """


rule cat_real_contigs:
    input:
        fa1=pjoin(WD, "input", "asm-{a}", "hap1.fa"),
        fa2=pjoin(WD, "input", "asm-{a}", "hap2.fa"),
    output:
        fa=pjoin(WD, "input", "asm-{a}", "haps.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule cat_alternative_contigs:
    input:
        fa1=pjoin(WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.hap1-w{w}.fa"),
        fa2=pjoin(WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.hap2-w{w}.fa"),
    output:
        fa=pjoin(WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.haps-w{w}.fa"),
    shell:
        """
        cat {input.fa1} {input.fa2} > {output.fa}
        """


rule align_alternative_contigs:
    input:
        tfa=pjoin(WD, "input", "asm-{a}", "haps.fa"),
        qfa=pjoin(WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.haps-w{w}.fa"),
    output:
        paf=pjoin(
            WD, "{ref}", "asmcallsets-{a}{confident}", "{asmcaller}.haps-w{w}.paf"
        ),
    conda:
        "../envs/minimap2.yml"
    threads: workflow.cores
    shell:
        """
        minimap2 -t{threads} -c {input.tfa} {input.qfa} > {output.paf}
        """
