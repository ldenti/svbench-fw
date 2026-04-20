import sys
import os
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

x_order = [
    "SVDSS",
    # "SVDSS2",
    "cuteSV",
    "debreak",
    "sawfish",
    "severus",
    "sniffles",
    "SVision-pro",
]

sns.set(style="whitegrid")

ref_map = {"chm13": "T2T-CHM13", "hg38": "GRCh38", "hg19": "GRCh37"}
cal_map = {
    "SVDSS": "SVDSS",
    "cutesv": "cuteSV",
    "debreak": "debreak",
    "sawfish": "sawfish",
    "severus": "severus",
    "sniffles": "sniffles",
    "svisionpro": "SVision-pro",
}


def clean_df(df, confident, refined, supp, mapq):
    df["Reference"] = df["Reference"].replace(ref_map)
    df["Caller"] = df["Caller"].replace(cal_map)

    conf = "strat-conf" if confident else "Full"

    df = df[(df["Setting"] == conf) & (df["Refined"] == refined)]

    sub_dfs = []
    for tool in df["Caller"].unique():
        if tool == "SVDSS":
            s, m = 4, 0
        elif tool == "cuteSV":
            s, m = 4, 20
        elif tool == "debreak":
            s, m = 0, 20
        elif tool == "sawfish":
            s, m = 0, 20
        elif tool == "severus":
            s, m = 4, 20
        elif tool == "sniffles":
            s, m = 0, 20
        elif tool == "SVision-pro":
            s, m = 4, 20
        else:
            assert False

        if supp != 0:
            if tool != "sawfish":
                s = supp
        if mapq != 0:
            if tool != "SVDSS":
                m = mapq

        sub_dfs.append(
            df[(df["Caller"] == tool) & (df["Support"] == s) & (df["MAPQ"] == m)]
        )
    return pd.concat(sub_dfs, axis=0, ignore_index=True)


# # Supplementary Table 5 (full table)
# df.sort_values(["RefSeq", "Truth", "Bench", "Refine"]).to_csv(
#     sys.stdout, index=False
# )

# # Supplementary Table 6 (delta(refine,norefine)
# print(
#     "RefSeq", "Truth", "Bench", "Tool", "F1-refine", "F1-norefine", "delta", sep=","
# )
# for refseq in df["RefSeq"].unique():
#     for truth in df["Truth"].unique():
#         if truth == "v0.6" and refseq != "GRCh37":
#             continue
#         for bench in df["Bench"].unique():
#             for tool in df["Tool"].unique():
#                 avg_f1_refine = df[
#                     (df["RefSeq"] == refseq)
#                     & (df["Truth"] == truth)
#                     & (df["Bench"] == bench)
#                     & (df["Tool"] == tool)
#                     & (df["Refine"] == True)
#                 ]["F1"].iloc[0]
#                 avg_f1_norefine = df[
#                     (df["RefSeq"] == refseq)
#                     & (df["Truth"] == truth)
#                     & (df["Bench"] == bench)
#                     & (df["Tool"] == tool)
#                     & (df["Refine"] == False)
#                 ]["F1"].iloc[0]
#                 print(
#                     refseq,
#                     truth,
#                     bench,
#                     tool,
#                     avg_f1_refine,
#                     avg_f1_norefine,
#                     round(avg_f1_refine - avg_f1_norefine, 2),
#                     sep=",",
#                 )

# # Correlation between 06 and 11
#         print(f"### {bench} /", "Refine" if refine else "NoRefine", "###")
#         old_f1 = df[
#             (df["Bench"] == bench)
#             & (df["Refine"] == refine)
#             & (df["Truth"] == "v0.6")
#         ].sort_values(["Tool"])[["Tool", "F1"]]
#         new_f1 = df[
#             (df["Bench"] == bench)
#             & (df["Refine"] == refine)
#             & (df["Truth"] == "v1.1")
#             & (df["RefSeq"] == "GRCh37")
#         ].sort_values(["Tool"])[["Tool", "F1"]]
#         print(old_f1)
#         print(new_f1)
#         corr = np.corrcoef(old_f1["F1"], new_f1["F1"])[0, 1]
#         print(corr)


def main():
    sns.set(font_scale=0.75)

    parser = argparse.ArgumentParser()
    parser.add_argument("CSV")
    parser.add_argument("--confident", action="store_true")
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("-m", "--mapq", type=int, default=0)
    parser.add_argument("-s", "--supp", type=int, default=0)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)
    df = clean_df(df, args.confident, args.refine, args.mapq, args.supp)

    df = df[["Reference", "GIAB", "Caller", "P", "R", "F1"]]
    df["rank"] = (
        df.groupby(["Reference", "GIAB"])["F1"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )
    print(df.to_csv(index=False))

    df["P"] = df["P"] * 100
    df["R"] = df["R"] * 100
    df["F1"] = df["F1"] * 100

    fig, axes = plt.subplots(2, 4, sharey="row", figsize=(10, 5))

    # GIAB v0.6
    sub_df = df[df["GIAB"] == 6]
    # print(sub_df)
    ref = sub_df["Reference"].unique()[0]
    col = 0

    sns.scatterplot(
        data=sub_df,
        x="P",
        y="R",
        hue="Caller",
        hue_order=x_order,
        ax=axes[0][col],
        legend=None,
    )
    sns.barplot(
        data=sub_df,
        x="Caller",
        y="F1",
        order=x_order,
        hue_order=x_order,
        ax=axes[1][col],
        hue="Caller",
    )

    axes[0][col].set_title(
        f"(a)\n" + ref + " (GIAB v0.6" + (")" if not args.confident else ", Tier 1)")
    )

    axes[1][col].set_ylim([0, 100])
    for i, container in enumerate(axes[1][col].containers):
        axes[1][col].bar_label(
            container, labels=sub_df[sub_df["Caller"] == x_order[i]]["rank"]
        )
    axes[1][col].tick_params(axis="x", labelrotation=90)

    # GIAB v5.0q
    ref = sub_df["Reference"].unique()[0]
    for col, ref in enumerate(
        # df["Reference"].unique(),
        ["GRCh37", "GRCh38", "T2T-CHM13"],
        1,
    ):
        sub_df = df[(df["Reference"] == ref) & (df["GIAB"] == 50)]

        sns.scatterplot(
            data=sub_df,
            x="P",
            y="R",
            hue="Caller",
            hue_order=x_order,
            ax=axes[0][col],
            legend=None,
        )

        sns.barplot(
            data=sub_df,
            x="Caller",
            y="F1",
            order=x_order,
            hue_order=x_order,
            ax=axes[1][col],
            hue="Caller",
            # alpha=0.75,
        )

        axes[0][col].set_title(
            f"({('abcd'[col])})\n"
            + ref
            + " (GIAB v5.0q"
            + (")" if not args.confident else ", w/ BED)")
        )
        axes[1][col].set_ylim([0, 100])

        for i, container in enumerate(axes[1][col].containers):
            axes[1][col].bar_label(
                container,
                labels=sub_df[(sub_df["Caller"] == x_order[i])]["rank"],
            )

        axes[1][col].tick_params(axis="x", labelrotation=90)

        if col != 0:
            axes[0][col].set_ylabel("")
            axes[1][col].set_ylabel("")

    xlim = 40  # 80 if args.bed else 40
    ylim = 30  # 40 if args.bed else 30
    for col in range(4):
        axes[0][col].set_xlim([xlim, 100])
        axes[0][col].set_ylim([ylim, 100])
    plt.tight_layout()

    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
