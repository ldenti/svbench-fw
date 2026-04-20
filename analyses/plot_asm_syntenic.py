import argparse
import os
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


CALLER_MAP = {"dipcall": "dipcall", "hapdiff": "hapdiff", "svim-asm": "SVIM-asm"}


def parse_vcf(vcf_fn, caller, tag):
    data = []
    for record in VariantFile(vcf_fn):
        l = len(record.alts[0]) - len(record.ref)
        gt1, gt2 = record.samples[0]["GT"]
        gt1 = gt1 if gt1 != None else 0
        gt2 = gt2 if gt2 != None else 0

        if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
            data.append(
                [
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
        df += parse_vcf(
            os.path.join(
                args.WD, "chm13", f"asmcallsets-{asm}{conf}", f"{caller}.vcf.gz"
            ),
            CALLER_MAP[caller],
            "Total",
        )
        df += parse_vcf(
            os.path.join(
                args.WD,
                "chm13",
                f"asmcallsets-{asm}{conf}",
                f"{caller}.non-syntenic.vcf",
            ),
            CALLER_MAP[caller],
            "Non-syntenic",
        )
    df = pd.DataFrame(df, columns=["Caller", "Length", "Type", "Tag"])
    # print(df)

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(9, 4))
    for col, caller in enumerate(["dipcall", "SVIM-asm", "hapdiff"]):
        ax = axes[col]
        subdf = df[(df["Caller"] == caller) & (abs(df["Length"]) <= 500)]

        sns.histplot(
            data=subdf,  # [subdf["Tag"] == "Syntenic"],
            x="Length",
            hue="Tag",
            # hue_order=TRUTHS,
            element="poly",
            fill=True,
            legend=True if col == 1 else None,
            bins=100,
            ax=ax,
        )
        ax.set_title(caller)
        if col == 0:
            ax.set_ylabel("Count")
        if col == 1:
            sns.move_legend(
                ax,
                "upper left",
                # bbox_to_anchor=(0.5, 1.02),
                ncol=1,
                title=None,  # "Regions",
                frameon=False,
                # columnspacing=0.1,
                # labelspacing=0.01,
            )
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
