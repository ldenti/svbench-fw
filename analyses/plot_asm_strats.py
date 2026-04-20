import argparse
import os
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

CALLER_MAP = {"dipcall": "dipcall", "hapdiff": "hapdiff", "svim-asm": "SVIM-asm"}
STRAT_MAP = {
    "easy": "Easy",
    "lowcomplex": "LowCpx",
    "lowmap": "LowMap",
    "segdups": "SegDups",
}


def parse_vcf(vcf_fn, asm, caller, tag):
    data = []
    for record in VariantFile(vcf_fn):
        l = len(record.alts[0]) - len(record.ref)
        gt1, gt2 = record.samples[0]["GT"]
        gt1 = gt1 if gt1 != None else 0
        gt2 = gt2 if gt2 != None else 0

        if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
            data.append(
                [
                    asm,
                    caller,
                    l,
                    "INS" if l > 0 else "DEL",
                    tag,
                    # gt1 == gt2,
                ]
            )
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", default="")
    parser.add_argument("--asm", required=True, type=str)
    parser.add_argument("--confident", action="store_true")
    parser.add_argument("WD")

    args = parser.parse_args()

    conf = ".confident" if args.confident else ""

    df = []
    asm = args.asm
    for caller in ["dipcall", "svim-asm", "hapdiff"]:
        # df += parse_vcf(
        #     os.path.join(
        #         args.WD, "chm13", f"asmcallsets-{asm}", f"{caller}.vcf.gz"
        #     ),
        #     asm,
        #     caller,
        #     "All",
        # )
        for strat in ["easy", "lowcomplex", "lowmap", "segdups"]:
            vcf_fn = os.path.join(
                args.WD, "chm13", f"asmcallsets-{asm}" + conf, f"{caller}.{strat}.vcf"
            )
            df += parse_vcf(vcf_fn, asm, CALLER_MAP[caller], STRAT_MAP[strat])
    df = pd.DataFrame(df, columns=["Assembly", "Caller", "Length", "Type", "Region"])
    print(df)

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(9, 4))
    for col, caller in enumerate(["dipcall", "SVIM-asm", "hapdiff"]):
        ax = axes[col]
        subdf = df[
            (df["Assembly"] == asm)
            & (df["Caller"] == caller)
            & (abs(df["Length"]) <= 500)
        ]

        sns.histplot(
            data=subdf,
            x="Length",
            hue="Region",
            element="poly",
            fill=False,
            legend=True if col == 1 else None,
            bins=100,
            ax=ax,
        )
        if col == 1:
            sns.move_legend(
                ax,
                "upper center",
                bbox_to_anchor=(0.5, 1.02),
                ncol=2,
                title="Regions",
                frameon=False,
                columnspacing=0.5,
                labelspacing=0.1,
            )
        ax.set_title(caller)
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
