rule truvari:
    input:
        fa=REF,
        vcf=pjoin(WD, "callsets", "{caller}.vcf.gz"),
        truth=pjoin(WD, "truths", "{truth}.vcf.gz"),
        bed=lambda wildcards: HAPBED if wildcards.truth == "hapdiff" else DIPBED,
    output:
        d=directory(pjoin(WD, "truvari-{truth}-{option}", "{caller}")),
        dd=directory(pjoin(WD, "truvari-{truth}-{option}", "{caller}", "phab_bench")),
    params:
        opt=lambda wildcards, input: truvari_options[wildcards.option]
        + (" " + input.bed if wildcards.option == "wbed" else ""),
        aligner="mafft",  # lambda wildcards: "mafft" if wildcards.option == "def" else "poa",
    threads: workflow.cores / 2
    conda:
        "../envs/truvari.yml"
    shell:
        """
        rm -rf {output.d}
        truvari bench {params.opt} --reference {input.fa} --base {input.truth} --comp {input.vcf} --output {output.d}
        truvari refine --reference {input.fa} --regions {output.d}/candidate.refine.bed --coords R --use-original-vcfs --threads {threads} --align {params.aligner} {output.d}
        # we need this if since this will fail if there are not regions to refine (like in the example data)
        # wasn't an issue with real data
        if [ -d {output.dd} ] ; then truvari ga4gh --input {output.d} --output {output.d}/ga4gh_with_refine ; else mkdir {output.dd} ; touch {output.dd}/NOREGIONSTOREFINE ; fi # with-refine
        """


rule format_truvari:
    input:
        expand(pjoin(WD, "truvari-{{truth}}-{{option}}", "{caller}"), caller=CALLERS),
    output:
        csv=pjoin(WD, "{truth}.truvari-{option}.csv"),
        refcsv=pjoin(WD, "{truth}.truvari-{option}.refine.csv"),
    params:
        bd=pjoin(WD, "truvari-{truth}-{option}"),
    shell:
        """
        python3 ./scripts/format_truvari.py {params.bd} > {output.csv}
        python3 ./scripts/format_truvari.py --refine {params.bd} > {output.refcsv}
        """


rule plot_truvari:
    input:
        pjoin(WD, "{truth}.truvari-{option}.csv"),
    output:
        pjoin(WD, "{truth}.truvari-{option}.csv.png"),
    conda:
        "../envs/seaborn.yml"
    shell:
        """
        python3 ./scripts/plot_truvari.py single {input}
        """
