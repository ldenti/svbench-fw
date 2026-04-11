# from snakemake.utils import min_version
# min_version("6.4.1")
import sys, os
from os.path import join as pjoin


##### config file #####
# configfile: "config/config.yml"


def lns(src, dst):
    src = os.path.realpath(src)
    if os.path.exists(dst):
        # print(os.readlink(dst), src)
        if os.readlink(dst) == src:
            return
        print(f"ln -s -f {src} {dst}")
        os.remove(dst)
    else:
        print(f"ln -s {src} {dst}")
    os.symlink(src, dst)


WD = config["wd"]

os.makedirs(WD, exist_ok=True)
os.makedirs(pjoin(WD, "input"), exist_ok=True)

# REFs = config["ref"]
references = []
os.makedirs(pjoin(WD, "input", "refs"), exist_ok=True)
for ref, fn in config["ref"].items():
    references.append(ref)
    # XXX: assuming not gzipped, but it should still work if no tool makes any assumption from extension
    lns(fn, pjoin(WD, "input", "refs", f"{ref}.fa"))
    if not os.path.exists(fn + ".fai"):
        print(
            "ERROR: please index the reference with 'samtools faidx'",
            file=sys.stderr,
        )
        sys.exit(1)
    lns(fn + ".fai", pjoin(WD, "input", "refs", f"{ref}.fa.fai"))

# TRFs = config["trf"]
os.makedirs(pjoin(WD, "input", "trfs"), exist_ok=True)
for ref, fn in config["trf"].items():
    lns(fn, pjoin(WD, "input", "trfs", f"{ref}.bed"))

# STRATs = config["strat"]
strats = []
os.makedirs(pjoin(WD, "input", "strats"), exist_ok=True)
for strat, fns in config["strat"].items():
    strats.append(strat)
    os.makedirs(pjoin(WD, "input", "strats", strat), exist_ok=True)
    for ref, fn in fns.items():
        lns(fn, pjoin(WD, "input", "strats", strat, f"{ref}.bed"))


if "giab11" in config:
    for t, fns in config["giab11"].items():
        os.makedirs(pjoin(WD, "input", "giab11"), exist_ok=True)
        for ref, fn in fns.items():
            if t == "bed":
                lns(fn, pjoin(WD, "input", "giab11", f"{ref}.bed"))
            else:
                lns(fn, pjoin(WD, "input", "giab11", f"{ref}.vcf.gz"))
                lns(fn + ".tbi", pjoin(WD, "input", "giab11", f"{ref}.vcf.gz.tbi"))

if "giab06" in config:
    for t, fns in config["giab06"].items():
        os.makedirs(pjoin(WD, "input", "giab06"), exist_ok=True)
        for ref, fn in fns.items():
            if t == "bed":
                lns(fn, pjoin(WD, "input", "giab06", f"{ref}.bed"))
            else:
                lns(fn, pjoin(WD, "input", "giab06", f"{ref}.vcf.gz"))
                lns(fn + ".tbi", pjoin(WD, "input", "giab06", f"{ref}.vcf.gz.tbi"))

SAMPLE_NAME = config["name"]
FQ = config["fq"]

asms = []
for n, fns in config["haps"].items():
    asms.append(n)
    os.makedirs(pjoin(WD, "input", f"asm-{n}"), exist_ok=True)
    for i, fn in enumerate(fns, 1):
        # XXX: assuming not gzipped, but it should still work if no tool makes any assumption from extension
        lns(fn, pjoin(WD, "input", f"asm-{n}", f"hap{i}.fa"))

ASMC = config["asm-callers"]
TAT_STRINGS = [f"{x}-against-{y}" for x in ASMC for y in ASMC if x != y]
print("===", len(ASMC), "ASSEMBLY-BASED CALLERS: ", ", ".join(ASMC))

CALLERS = []
for caller, opts in config["callers"].items():
    for s in opts["s"]:
        for q in opts["q"]:
            CALLERS.append(f"{caller}-s{s}-q{q}")
print("===", len(CALLERS), "READ-BASED CALLERS:", ", ".join(CALLERS))

"""
A note on callers: I tried to use output VCFs as they are. But sometimes,
truvari was complaining about BND records, so I decided to remove them.
"""


wildcard_constraints:
    ref="|".join([f"({x})" for x in references]),
    a="|".join([f"({x})" for x in asms]),
    s=r"\d+",
    q=r"\d+",
    strat="|".join([f"({x})" for x in strats]),
    asmcaller="|".join([f"({c})" for c in ASMC]),


# callers from assembly
include: "rules/callers-asm.smk"
include: "rules/callers-asm-post.smk"
# alignment and phasing
include: "rules/preprocess.smk"
# callers from BAM
include: "rules/cutesv.smk"
include: "rules/debreak.smk"
include: "rules/sawfish.smk"
include: "rules/severus.smk"
include: "rules/sniffles.smk"
include: "rules/svision-pro.smk"
include: "rules/svdss2-ht.smk"
include: "rules/svdss.smk"
# truvari
include: "rules/truvari.smk"


rule all:
    input:
        # from callers-asm.smk
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmc}.vcf.gz"),
            ref=references,
            a=asms,
            asmc=ASMC,
        ),
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}.confident", "{asmc}.vcf.gz"),
            ref=references,
            a=asms,
            asmc=ASMC,
        ),
        # from various callers .smk
        expand(
            pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
            ref=references,
            caller=CALLERS,
        ),
        #
        pjoin(WD, "truvari.csv"),
        pjoin(WD, "truvari-giab.csv"),
        #
        expand(
            pjoin(
                WD,
                "{ref}",
                "asmcallsets-{a}",
                "comparison-{mode}",
                "{tat}",
            ),
            ref=references,
            a=asms,
            mode=["full", "conf"],
            tat=TAT_STRINGS,
        ),
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmc}.haps-w{w}.paf"),
            ref=references,
            a=asms,
            asmc=ASMC,
            w=[500],
        ),
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}.confident", "{asmc}.haps-w{w}.paf"),
            ref=references,
            a=asms,
            asmc=ASMC,
            w=[500],
        ),


rule asmcallers:
    input:
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmc}.vcf.gz"),
            ref=references,
            a=asms,
            asmc=ASMC,
        ),
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}", "{asmc}.haps-w{w}.paf"),
            ref=references,
            a=asms,
            asmc=ASMC,
            w=[500],
        ),
        expand(
            pjoin(WD, "{ref}", "asmcallsets-{a}.confident", "{asmc}.haps-w{w}.paf"),
            ref=references,
            a=asms,
            asmc=ASMC,
            w=[500],
        ),


rule asmcallers_comparison:
    input:
        expand(
            pjoin(
                WD,
                "{ref}",
                "asmcallsets-{a}",
                "comparison-{mode}",
                "{tat}",
            ),
            ref=references,
            a=asms,
            mode=["full", "conf"],
            tat=TAT_STRINGS,
        ),


rule callers:
    input:
        expand(
            pjoin(WD, "{ref}", "callsets", "{caller}.vcf.gz"),
            ref=references,
            caller=CALLERS,
        ),


rule asm_benchmarking_nostrats:
    input:
        # from truvari.smk
        expand(
            pjoin(WD, "{ref}", "truvari", "{x}", "{a}", "{asmc}", "{caller}"),
            ref=references,
            x=["full", "conf"],
            a=asms,
            asmc=ASMC,
            caller=CALLERS,
        ),
    output:
        pjoin(WD, "truvari.csv"),
    shell:
        """
        python3 ./scripts/format_truvari.py {WD} > {output}
        """


rule asm_benchmarking:
    input:
        # from truvari.smk
        expand(
            pjoin(WD, "{ref}", "truvari", "{x}", "{a}", "{asmc}", "{caller}"),
            ref=references,
            x=["full", "conf"] + strats,
            a=asms,
            asmc=ASMC,
            caller=CALLERS,
        ),
    output:
        pjoin(WD, "truvari.csv"),
    shell:
        """
        python3 ./scripts/format_truvari.py {WD} > {output}
        """


rule giab_benchmarking:
    input:
        # from truvari.smk
        expand(
            pjoin(WD, "{ref}", "truvari-giab", "11", "{opt}", "{caller}"),
            ref=references,
            opt=["full", "conf"],
            caller=CALLERS,
        ),
        expand(
            pjoin(WD, "hg19", "truvari-giab", "06", "{opt}", "{caller}"),
            opt=["full", "conf"],
            caller=CALLERS,
        ),
    output:
        pjoin(WD, "truvari-giab.csv"),
    shell:
        """
        python3 ./scripts/format_truvari_giab.py {WD} > {output}
        """
