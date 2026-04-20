import sys
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# sns.set_theme(style="whitegrid")

asm_map = {"dipcall": "dipcall", "hapdiff": "hapdiff", "svim-asm": "SVIM-asm"}
bench_map = {
    # "hard": "Hard",
    "strat-lowcomplex": "LowComplexity",
    "strat-easy": "Easy",
    "strat-segdups": "SegDups",
    "strat-lowmap": "LowMappability",
}
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


def clean_df(df, asm, refined, mapq, supp):
    if asm in df["Assembly"].values:
        df = df[df["Assembly"] == asm]
    else:
        print(f"{asm} error", file=sys.stderr)

    df = df[df["Refined"] == refined]
    df = df[df["Setting"].isin(bench_map)]

    df["Reference"] = df["Reference"].replace(ref_map)
    df["ASM-Caller"] = df["ASM-Caller"].replace(asm_map)
    df["Caller"] = df["Caller"].replace(cal_map)
    df["Setting"] = df["Setting"].replace(bench_map)

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


def main():
    sns.set_theme(font_scale=0.8)

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--asm", type=str, default="hprc")
    parser.add_argument("--refine", action="store_true")
    parser.add_argument("-m", "--mapq", type=int, default=0)
    parser.add_argument("-s", "--supp", type=int, default=0)
    parser.add_argument("-o", type=str, default="")
    parser.add_argument("CSV")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)

    df = clean_df(df, args.asm, args.refine, args.mapq, args.supp)
    df = df[["Reference", "Setting", "ASM-Caller", "Caller", "P", "R", "F1"]]

    print(df.to_csv(index=False))

    df["F1"] = df["F1"] * 100

    tools = df["Caller"].unique()
    tools.sort()
    xticks = {tool: i for i, tool in enumerate(tools, 1)}

    # df = df[df["Refine"] == args.refine]

    xoff = 0.25
    x_offsets = {"dipcall": (-xoff, "r"), "SVIM-asm": (0, "g"), "hapdiff": (xoff, "b")}

    markers = {
        # "Hard": "x",
        "Easy": "o",
        "LowComplexity": "^",
        "LowMappability": "s",
        "SegDups": "x",
    }

    refseqs = ["GRCh37", "GRCh38", "T2T-CHM13"]

    fig, axes = plt.subplots(1, 3, figsize=(11, 5), sharey=True)

    for i, refseq in enumerate(refseqs):
        xticksvalues = []
        xtickslabels = []
        for tool, x in xticks.items():
            for truth, (offset, color) in x_offsets.items():
                for bench, marker in markers.items():
                    f1 = df[
                        (df["Reference"] == refseq)
                        & (df["ASM-Caller"] == truth)
                        & (df["Setting"] == bench)
                        & (df["Caller"] == tool)
                    ]["F1"].iloc[0]
                    # print(tool, x, truth, offset, bench, f1)
                    x1 = x + offset
                    axes[i].plot(x1, f1, f"{color}{marker}")

                if offset == 0:
                    xticksvalues.append(x)
                    xtickslabels.append(tool)

        axes[i].set_xticks(xticksvalues, xtickslabels, rotation=30)
        axes[i].set_title(refseq)
        axes[i].set_ylim(10, 103)
        if i == 0:
            axes[i].set_ylabel("F1")

    legend_elements = [
        Patch(facecolor="r", edgecolor="r", label="dipcall"),
        Patch(facecolor="g", edgecolor="g", label="SVIM-asm"),
        Patch(facecolor="b", edgecolor="b", label="hapdiff"),
    ]
    legend_elements += [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="black",
            label=label,
            linewidth=0,
            markerfacecolor="black",
            # markersize=15,
        )
        for label, marker in markers.items()
    ]
    axes[0].legend(handles=legend_elements)

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
