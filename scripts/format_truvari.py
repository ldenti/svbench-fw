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

    print(
        "Reference",
        "Assembly",
        "ASM-Caller",
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
        os.path.join(indir, "*", "truvari", "*", "*", "*", "summary.json")
    ):
        fields = json.split("/")

        tool = fields[-2]
        caller, s, q = tool.split("-")
        s = s[1:]
        q = q[1:]
        asmc = fields[-3]
        asm = fields[-4]
        ref = fields[-6]

        base_d = os.path.dirname(json)

        # FULL - bench
        tp_base = count_vcf(base_d + "/tp-base.vcf.gz")
        tp_comp = count_vcf(base_d + "/tp-comp.vcf.gz")
        fp = count_vcf(base_d + "/fp.vcf.gz")
        fn = count_vcf(base_d + "/fn.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            asm,
            asmc,
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

        # FULL - refine
        tp_base, fn = parse_vcf(base_d + "/ga4gh_with_refine.base.vcf.gz")
        tp_comp, fp = parse_vcf(base_d + "/ga4gh_with_refine.comp.vcf.gz")
        P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
        print(
            ref,
            asm,
            asmc,
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

        for vcf in glob.glob(os.path.join(base_d, "*", "tp-base.vcf.gz")):
            strat = vcf.split("/")[-2]
            if strat == "phab_bench":
                continue

            strat_d = os.path.dirname(vcf)

            # STRAT - bench
            tp_base = count_vcf(strat_d + "/tp-base.vcf.gz")
            tp_comp = count_vcf(strat_d + "/tp-comp.vcf.gz")
            fp = count_vcf(strat_d + "/fp.vcf.gz")
            fn = count_vcf(strat_d + "/fn.vcf.gz")
            P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
            print(
                ref,
                asm,
                asmc,
                strat,
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

            # STRAT - refine
            tp_base, fn = parse_vcf(strat_d + "/ga4gh_with_refine.base.vcf.gz")
            tp_comp, fp = parse_vcf(strat_d + "/ga4gh_with_refine.comp.vcf.gz")
            P, R, F = compute_prf(tp_comp, tp_base, fp, fn)
            print(
                ref,
                asm,
                asmc,
                strat,
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
