import sys
import os
import glob


def parse_summary(fpath):
    TP, FP, FN = 0, 0, 0
    for line in open(fpath):
        line = line.strip('\n \t"')
        if line.startswith('TP-base"'):
            if TP == 0:
                TP = int(line.split(" ")[1][:-1])
        elif line.startswith('FP"'):
            if FP == 0:
                FP = int(line.split(" ")[1][:-1])
        elif line.startswith('FN"'):
            if FN == 0:
                FN = int(line.split(" ")[1][:-1])
    return TP, FP, FN


def main():
    indir = sys.argv[1]

    print(
        "Reference",
        "GIAB",
        "Setting",
        "Refined",
        "Caller",
        "Support",
        "MAPQ",
        "TP",
        "FP",
        "FN",
        "P",
        "R",
        "F1",
        sep=",",
    )

    for fn in glob.glob(
        os.path.join(indir, "*", "truvari-giab", "*", "*", "*", "summary.json")
    ):
        fields = fn.split("/")

        tool = fields[-2]
        caller, s, q = tool.split("-")
        s = s[1:]
        q = q[1:]
        setting = fields[-3]
        giabv = fields[-4]
        ref = fields[-6]

        for refine in ["norefined", "refined"]:
            f = fn
            if refine == "refined":
                f = "/".join(fields[:-1]) + "/" + "ga4gh_with_refine.summary.json"
            tp, fp, fn = parse_summary(f)
            P = round(tp / (tp + fp) * 100, 1) if tp + fp > 0 else 0.0
            R = round(tp / (tp + fn) * 100, 1) if tp + fn > 0 else 0.0
            F = round(2 * (P * R) / (P + R) if P + R != 0 else 0.0, 1)
            print(
                ref,
                giabv,
                setting,
                refine == "refined",
                caller,
                s,
                q,
                tp,
                fp,
                fn,
                P,
                R,
                F,
                sep=",",
            )


if __name__ == "__main__":
    main()
