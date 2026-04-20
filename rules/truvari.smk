# options: --passonly --pick ac --dup-to-ins


rule truvari:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        asmcaller=pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.vcf.gz"),
    output:
        pjoin(
            WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "tp-base.vcf.gz"
        ),
        pjoin(
            WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "tp-comp.vcf.gz"
        ),
        pjoin(WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "fp.vcf.gz"),
        pjoin(WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "fn.vcf.gz"),
        pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
    params:
        wd=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
        ),
    threads: workflow.cores / 8
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari bench --passonly --pick ac --dup-to-ins --refine --reference {input.fa} --base {input.asmcaller} --comp {input.vcf} --output {params.wd}
        truvari ga4gh --input {params.wd} --output {params.wd}/ga4gh_with_refine
        """


# # we need this if since this will fail if there are not regions to refine (like in the example data)
# # wasn't an issue with real data
# if [ -d {output.dd} ]
# then
# truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine
# else
# mkdir -p {output.dd}
# touch {output.d}/NOREGIONSTOREFINE
# cp {output.d}/summary.json {output.d}/ga4gh_with_refine.summary.json
# fi


rule truvari_giab:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        vcf=pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
        truth=pjoin(WD, "input", "giab{giabv}", "{ref}.vcf.gz"),
    output:
        pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
    params:
        wd=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
        ),
    conda:
        "../envs/truvari.yml"
    threads: workflow.cores / 8
    shell:
        """
        truvari bench --passonly --pick ac --dup-to-ins --refine --reference {input.fa} --base {input.asmcaller} --comp {input.vcf} --output {params.wd}
        truvari ga4gh --input {params.wd} --output {params.wd}/ga4gh_with_refine
        """


rule stratify_giab:
    input:
        tp_base=pjoin(
            WD, "{ref}", "truvari-giab", "{giabv}", "{caller}", "tp-base.vcf.gz"
        ),
        tp_comp=pjoin(
            WD, "{ref}", "truvari-giab", "{giabv}", "{caller}", "tp-comp.vcf.gz"
        ),
        fp=pjoin(WD, "{ref}", "truvari-giab", "{giabv}", "{caller}", "fp.vcf.gz"),
        fn=pjoin(WD, "{ref}", "truvari-giab", "{giabv}", "{caller}", "fn.vcf.gz"),
        base=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        comp=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
        bed=pjoin(WD, "input", "giab{giabv}", "{ref}.bed"),
    output:
        tp_base=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "tp-base.vcf.gz",
        ),
        tp_comp=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "tp-comp.vcf.gz",
        ),
        fp=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "fp.vcf.gz",
        ),
        fn=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "fn.vcf.gz",
        ),
        base=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        comp=pjoin(
            WD,
            "{ref}",
            "truvari-giab",
            "{giabv}",
            "{caller}",
            "strat-conf",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -header -a {input.tp_base} -b {input.bed} -u | bgzip -c > {output.tp_base}
        tabix -p vcf {output.tp_base}
        bedtools intersect -header -a {input.tp_comp} -b {input.bed} -u | bgzip -c > {output.tp_comp}
        tabix -p vcf {output.tp_comp}
        bedtools intersect -header -a {input.fp} -b {input.bed} -u | bgzip -c > {output.fp}
        tabix -p vcf {output.fp}
        bedtools intersect -header -a {input.fn} -b {input.bed} -u | bgzip -c > {output.fn}
        tabix -p vcf {output.fn}
        #
        bedtools intersect -header -a {input.base} -b {input.bed} -u | bgzip -c > {output.base}
        tabix -p vcf {output.base}
        bedtools intersect -header -a {input.comp} -b {input.bed} -u | bgzip -c > {output.comp}
        tabix -p vcf {output.comp}
        """


rule stratify:
    input:
        tp_base=pjoin(
            WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "tp-base.vcf.gz"
        ),
        tp_comp=pjoin(
            WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "tp-comp.vcf.gz"
        ),
        fp=pjoin(WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "fp.vcf.gz"),
        fn=pjoin(WD, "{ref}", "truvari", "{a}", "{asmcaller}", "{caller}", "fn.vcf.gz"),
        base=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        comp=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
        bed=lambda wildcards: (
            pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmcaller}.bed")
            if wildcards.strat == "conf"
            else pjoin(WD, "input", "strats", "{strat}", "{ref}.bed")
        ),
    output:
        tp_base=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "tp-base.vcf.gz",
        ),
        tp_comp=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "tp-comp.vcf.gz",
        ),
        fp=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "fp.vcf.gz",
        ),
        fn=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "fn.vcf.gz",
        ),
        base=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "ga4gh_with_refine.base.vcf.gz",
        ),
        comp=pjoin(
            WD,
            "{ref}",
            "truvari",
            "{a}",
            "{asmcaller}",
            "{caller}",
            "strat-{strat}",
            "ga4gh_with_refine.comp.vcf.gz",
        ),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -header -a {input.tp_base} -b {input.bed} -u | bgzip -c > {output.tp_base}
        tabix -p vcf {output.tp_base}
        bedtools intersect -header -a {input.tp_comp} -b {input.bed} -u | bgzip -c > {output.tp_comp}
        tabix -p vcf {output.tp_comp}
        bedtools intersect -header -a {input.fp} -b {input.bed} -u | bgzip -c > {output.fp}
        tabix -p vcf {output.fp}
        bedtools intersect -header -a {input.fn} -b {input.bed} -u | bgzip -c > {output.fn}
        tabix -p vcf {output.fn}
        #
        bedtools intersect -header -a {input.base} -b {input.bed} -u | bgzip -c > {output.base}
        tabix -p vcf {output.base}
        bedtools intersect -header -a {input.comp} -b {input.bed} -u | bgzip -c > {output.comp}
        tabix -p vcf {output.comp}
        """
