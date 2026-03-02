# DIPBED = pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.bed"),
# HAPBED = pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff", "confident_regions.bed"),

# options: --passonly --pick ac --dup-to-ins


rule truvari_full:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        asmcaller=pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.vcf.gz"),
    output:
        d=directory(
            pjoin(WD, "{ref}", "truvari", "full", "{a}", "{asmcaller}", "{caller}")
        ),
        dd=directory(
            pjoin(
                WD,
                "{ref}",
                "truvari",
                "full",
                "{a}",
                "{asmcaller}",
                "{caller}",
                "phab_bench",
            )
        ),
    threads: workflow.cores / 4
    conda:
        "../envs/truvari.yml"
    shell:
        """
        # set +e
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins --reference {input.fa} --base {input.asmcaller} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        # we need this if since this will fail if there are not regions to refine (like in the example data)
        # wasn't an issue with real data
        if [ -d {output.dd} ]
        then
            truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine
        else
            mkdir -p {output.dd}
            touch {output.d}/NOREGIONSTOREFINE
            cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
        fi
        """


rule truvari_confident:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        asmcaller=pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.vcf.gz"),
        bed=lambda wildcards: (
            pjoin(WD, "{ref}", "asmcallsets-{a}", "hapdiff", "confident_regions.bed")
            if wildcards.asmcaller == "hapdiff"
            else pjoin(WD, "{ref}", "asmcallsets-{a}", "dipcall.bed")
        ),
    output:
        d=directory(
            pjoin(WD, "{ref}", "truvari", "conf", "{a}", "{asmcaller}", "{caller}")
        ),
        dd=directory(
            pjoin(
                WD,
                "{ref}",
                "truvari",
                "conf",
                "{a}",
                "{asmcaller}",
                "{caller}",
                "phab_bench",
            )
        ),
    threads: workflow.cores / 4
    conda:
        "../envs/truvari.yml"
    shell:
        """
        # set +e
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins --includebed {input.bed} --reference {input.fa} --base {input.asmcaller} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        # we need this if since this will fail if there are not regions to refine (like in the example data)
        # wasn't an issue with real data
        if [ -d {output.dd} ]
        then
            truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine
        else
            mkdir -p {output.dd}
            touch {output.d}/NOREGIONSTOREFINE
            cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
        fi
        """


rule truvari_strat:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        asmcaller=pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.vcf.gz"),
        bed=pjoin(WD, "input", "strats", "{strat}", "{ref}.bed"),
    output:
        d=directory(
            pjoin(WD, "{ref}", "truvari", "{strat}", "{a}", "{asmcaller}", "{caller}")
        ),
        dd=directory(
            pjoin(
                WD,
                "{ref}",
                "truvari",
                "{strat}",
                "{a}",
                "{asmcaller}",
                "{caller}",
                "phab_bench",
            )
        ),
    threads: workflow.cores / 4
    conda:
        "../envs/truvari.yml"
    shell:
        """
        # set +e
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins --includebed {input.bed} --reference {input.fa} --base {input.asmcaller} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        # we need this if since this will fail if there are not regions to refine (like in the example data)
        # wasn't an issue with real data
        if [ -d {output.dd} ]
        then
            truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine
        else
            mkdir -p {output.dd}
            touch {output.d}/NOREGIONSTOREFINE
            cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
        fi
        """

rule truvari_giab:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        truth=pjoin(WD, "input", "giab{giabv}", "{ref}.vcf.gz"),
        bed=pjoin(WD, "input", "giab{giabv}", "{ref}.bed"),
    output:
        d=directory(pjoin(WD, "{ref}", "truvari-giab", "{giabv}", "{opt}", "{caller}")),
        dd=directory(
            pjoin(
                WD,
                "{ref}",
                "truvari-giab",
                "{giabv}",
                "{opt}",
                "{caller}",
                "phab_bench",
            )
        ),
    params:
        bed=lambda wildcards, input: (
            "--includebed " + input.bed if wildcards.opt == "conf" else ""
        ),
    conda:
        "../envs/truvari.yml"
    threads: workflow.cores / 4
    shell:
        """
        # set +e
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins {params.bed} --reference {input.fa} --base {input.truth} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        # we need this if since this will fail if there are not regions to refine (like in the example data)
        # wasn't an issue with real data
        if [ -d {output.dd} ]
        then
            truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine # --with-refine
        else
            mkdir -p {output.dd}
            touch {output.d}/NOREGIONSTOREFINE
            cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
        fi
        """
