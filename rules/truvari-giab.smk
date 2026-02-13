# wildcard_constraints:
#     giab="(11)|(06)",
#     giabopt="(full)|(conf)",


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
    threads: workflow.cores / 2
    shell:
        """
        set +e
        rm -rf {output.d}
        truvari bench --passonly --pick ac --dup-to-ins {params.bed} --reference {input.fa} --base {input.truth} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align mafft {output.d}
        exitcode=$?
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
