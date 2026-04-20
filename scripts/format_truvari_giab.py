import sys
import os
import glob

from pysam import VariantFile


def count_vcf(fn):
    n = 0
    for record in VariantFile(fn):
        n += 1
    return n


def parse_vcf(fn):
    tp, other = 0, 0
    for record in VariantFile(fn):
        bd = record.samples[0]["BD"]
        if bd == "TP":
            tp += 1
        else:
            assert bd in ["FN", "FP"]
            other += 1
    return tp, other


def compute_prf(tp_comp, tp_base, fp, fn):
    P = tp_comp / (tp_comp + fp)
    R = tp_base / (tp_base + fn)
    F = 2 * P * R / (P + R)
    return P, R, F


def main():
    indir = sys.argv[1]

    CALLERS = [
        "debreak-s0-q20",
        "sawfish-s0-q20",
        "sniffles-s0-q20",
        "cutesv-s4-q20",
        "severus-s4-q20",
        "svisionpro-s4-q20",
        "SVDSS-s4-q0",
    ]

    print(
        "Reference",
        "GIAB",
        "Setting",
        "Refined",
        "Caller",
        "Support",
        "MAPQ",
        "TP-base",
        "TP-comp",
        "FP",
        "FN",
        "P",
        "R",
        "F1",
        sep=",",
    )

    for json in glob.glob(
        os.path.join(indir, "*", "truvari-giab", "*", "*", "summary.json")
    ):
        fields = json.split("/")

        tool = fields[-2]
        if tool not in CALLERS:
            continue
        caller, s, q = tool.split("-")
        s = s[1:]
        q = q[1:]
        giabv = fields[-3]
        ref = fields[-5]

        base_d = os.path.dirname(json)

        # FULL GENOME (bench)
        tp_base = count_vcf(base_d + "/tp-base.vcf.gz")
        tp_comp = count_vcf(base_d + "/tp-comp.vcf.gz")
        fp = count_vcf(base_d + "/fp.vcf.gz")
        fn = count_vcf(base_d + "/fn.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            giabv,
            "Full",
            "False",
            caller,
            s,
            q,
            tp_base,
            tp_comp,
            fp,
            fn,
            P,
            R,
            F,
            sep=",",
        )

        # FULL GENOME (refine)
        tp_base, fn = parse_vcf(base_d + "/ga4gh_with_refine.base.vcf.gz")
        tp_comp, fp = parse_vcf(base_d + "/ga4gh_with_refine.comp.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            giabv,
            "Full",
            "True",
            caller,
            s,
            q,
            tp_base,
            tp_comp,
            fp,
            fn,
            P,
            R,
            F,
            sep=",",
        )

        strat_d = os.path.join(base_d, "strat-conf")

        # CONFIDENT REGIONS (bench)
        tp_base = count_vcf(strat_d + "/tp-base.vcf.gz")
        tp_comp = count_vcf(strat_d + "/tp-comp.vcf.gz")
        fp = count_vcf(strat_d + "/fp.vcf.gz")
        fn = count_vcf(strat_d + "/fn.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            giabv,
            "strat-conf",
            "False",
            caller,
            s,
            q,
            tp_base,
            tp_comp,
            fp,
            fn,
            P,
            R,
            F,
            sep=",",
        )

        # CONFIDENT REGIONS (refine)
        tp_base, fn = parse_vcf(strat_d + "/ga4gh_with_refine.base.vcf.gz")
        tp_comp, fp = parse_vcf(strat_d + "/ga4gh_with_refine.comp.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            giabv,
            "strat-conf",
            "True",
            caller,
            s,
            q,
            tp_base,
            tp_comp,
            fp,
            fn,
            P,
            R,
            F,
            sep=",",
        )


if __name__ == "__main__":
    main()
