import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="whitegrid")

ref_map = {"chm13": "T2T-CHM13", "hg38": "GRCh38", "hg19": "GRCh37"}
asm_map = {"dipcall": "dipcall", "hapdiff": "hapdiff", "svim-asm": "SVIM-asm"}
cal_map = {
    "SVDSS": "SVDSS",
    "cutesv": "cuteSV",
    "debreak": "debreak",
    "sawfish": "sawfish",
    "severus": "severus",
    "sniffles": "sniffles",
    "svisionpro": "SVision-pro",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("CSV")
    parser.add_argument("-a", "--asm", required=True)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)

    df = df[df["Assembly"] == args.asm]
    df = df[df["Refined"] == True]

    df["Reference"] = df["Reference"].replace(ref_map)
    df["ASM-Caller"] = df["ASM-Caller"].replace(asm_map)
    df["Caller"] = df["Caller"].replace(cal_map)

    df = df[df["Setting"].isin(["Full", "strat-conf"])]

    # Supplementary Table
    for col, ref in enumerate(["GRCh37", "GRCh38", "T2T-CHM13"]):
        for row, region in enumerate(["Full", "strat-conf"]):
            subdf = df[(df["Reference"] == ref) & (df["Setting"] == region)]

            r = (
                subdf.sort_values(
                    ["ASM-Caller", "Caller", "F1"], ascending=[True, True, False]
                )
                .groupby(["ASM-Caller", "Caller"], sort=False)
                .head(1)
                .reset_index(drop=True)
            )
            for truth in subdf["ASM-Caller"].unique():
                for caller in subdf["Caller"].unique():
                    best = r[(r["ASM-Caller"] == truth) & (r["Caller"] == caller)][
                        "F1"
                    ].iloc[0]

                    s = 4
                    mq = 20
                    if caller == "SVDSS":
                        s, mq = 4, 0
                    elif caller == "cuteSV":
                        s, mq = 4, 20
                    elif caller == "debreak":
                        s, mq = 0, 20
                    elif caller == "sawfish":
                        s, mq = 0, 20
                    elif caller == "severus":
                        s, mq = 4, 20
                    elif caller == "sniffles":
                        s, mq = 0, 20
                    elif caller == "SVision-pro":
                        s, mq = 4, 20

                    selected = subdf[
                        (subdf["ASM-Caller"] == truth)
                        & (subdf["Caller"] == caller)
                        & (subdf["Support"] == s)
                        & (subdf["MAPQ"] == mq)
                    ]["F1"].iloc[0]

                    print(
                        ref,
                        region,
                        truth,
                        caller,
                        best,
                        selected,
                        best - selected,
                        sep=",",
                    )

    # Supplementary Figure
    df["F1"] = df["F1"] * 100

    fig, axes = plt.subplots(2, 3, figsize=(7, 7), sharex=True, sharey=True)

    for col, ref in enumerate(["GRCh37", "GRCh38", "T2T-CHM13"]):
        for row, region in enumerate(["Full", "strat-conf"]):
            subdf = df[(df["Reference"] == ref) & (df["Setting"] == region)]

            ax = axes[row][col]

            sns.boxplot(
                data=subdf,
                x="Caller",
                y="F1",
                hue="ASM-Caller",
                legend=True if row == 1 and col == 0 else False,
                ax=ax,
            )
            axes[row][col].set_ylim(40, 105)
        axes[0][col].set_title(ref)

        axes[1][col].set_xticklabels(axes[1][col].get_xticklabels(), rotation=90)
    axes[0][0].set_ylabel("Full Genome\nF1")
    axes[1][0].set_ylabel("Confident Regions\nF1")

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
