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


def clean_df(df, ref, asm, supp, mapq):
    if ref in df["Reference"].unique():
        df = df[df["Reference"] == ref]

    if asm in df["Assembly"].values:
        df = df[df["Assembly"] == asm]
    else:
        print(f"{asm} error", file=sys.stderr)

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


# Supplementary Table (delta(refine,norefine)
def supp_table():
    parser = argparse.ArgumentParser()
    parser.add_argument("CSV")
    parser.add_argument("-a", "--asm", type=str, default="hprc")
    parser.add_argument("-f", "--full", action="store_true")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)

    df["Reference"] = df["Reference"].replace(ref_map)
    df["ASM-Caller"] = df["ASM-Caller"].replace(asm_map)
    df["Caller"] = df["Caller"].replace(cal_map)

    if not args.full:
        df = clean_df(df, "all", ".", 0, 0)

    print(
        "Reference",
        "Assembly",
        "ASM-Caller",
        "Setting",
        "Caller",
        "Support",
        "MAPQ",
        "F1-refine",
        "F1-norefine",
        "delta",
        sep=",",
    )
    for ref in df["Reference"].unique():
        for asm in df["Assembly"].unique():
            for asmc in df["ASM-Caller"].unique():
                for bench in df["Setting"].unique():
                    for tool in df["Caller"].unique():
                        sub_df = df[
                            (df["Reference"] == ref)
                            & (df["Assembly"] == asm)
                            & (df["ASM-Caller"] == asmc)
                            & (df["Setting"] == bench)
                            & (df["Caller"] == tool)
                        ]
                        for supp in sub_df["Support"].unique():
                            for mapq in sub_df["MAPQ"].unique():
                                avg_f1_refine = sub_df[
                                    (sub_df["Support"] == supp)
                                    & (sub_df["MAPQ"] == mapq)
                                    & (sub_df["Refined"] == True)
                                ]["F1"].iloc[0]
                                avg_f1_norefine = sub_df[
                                    (sub_df["Support"] == supp)
                                    & (sub_df["MAPQ"] == mapq)
                                    & (sub_df["Refined"] == False)
                                ]["F1"].iloc[0]

                                print(
                                    ref,
                                    asm,
                                    asmc,
                                    bench,
                                    tool,
                                    supp,
                                    mapq,
                                    avg_f1_refine,
                                    avg_f1_norefine,
                                    round(avg_f1_refine - avg_f1_norefine, 2),
                                    sep=",",
                                )


def main_strip():
    parser = argparse.ArgumentParser()
    parser.add_argument("CSV")
    parser.add_argument("-a", "--asm", required=True)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("-m", "--mapq", type=int, default=0)
    parser.add_argument("-s", "--supp", type=int, default=0)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)
    print(df)

    df["Reference"] = df["Reference"].replace(ref_map)
    df["ASM-Caller"] = df["ASM-Caller"].replace(asm_map)
    df["Caller"] = df["Caller"].replace(cal_map)

    df = clean_df(df, "all", args.asm, args.supp, args.mapq)
    if not args.all:
        df = df[df["Setting"].isin(["Full", "strat-conf"])]

    # sns.set(font_scale=0.7)

    df["F1"] = df["F1"] * 100

    # Numbers for manuscript: avg f1
    for refseq in df["Reference"].unique():
        for bench in df["Setting"].unique():
            for refine in df["Refined"].unique():
                avg_f1 = df[
                    (df["Reference"] == refseq)
                    & (df["Setting"] == bench)
                    & (df["Refined"] == refine)
                ]["F1"].mean()
                print(refseq, bench, "ref" if refine else "noref", avg_f1, sep="\t")

    tools = list(df["Caller"].unique())
    tools.sort()

    print(df)

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6), sharey=True)
    sns.stripplot(
        df,
        x="F1",
        y="Caller",
        hue="Reference",
        alpha=0.6,
        s=6,
        order=tools,
        hue_order=ref_map.values(),
        linewidth=1,
        ax=ax1,
    )
    ax1.set_xlim(40, 100)
    sns.move_legend(
        ax1,
        "lower center",
        bbox_to_anchor=(0.5, 1),
        ncol=3,
        title=None,
        frameon=False,
    )

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


def main_heat():
    parser = argparse.ArgumentParser()
    parser.add_argument("CSV")
    parser.add_argument("-r", "--ref", type=str, default="chm13")
    parser.add_argument("-a", "--asm", type=str, default="hprc")
    parser.add_argument("-m", "--mapq", type=int, default=0)
    parser.add_argument("-s", "--supp", type=int, default=0)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()

    df = pd.read_csv(args.CSV)

    assert args.ref in ref_map
    ref_name = ref_map[args.ref]

    df["Reference"] = df["Reference"].replace(ref_map)
    df["ASM-Caller"] = df["ASM-Caller"].replace(asm_map)
    df["Caller"] = df["Caller"].replace(cal_map)

    df = clean_df(df, ref_name, args.asm, args.supp, args.mapq)

    truths = list(df["ASM-Caller"].unique())
    truths.sort()
    truths = ["dipcall", "SVIM-asm", "hapdiff"]

    tools = list(df["Caller"].unique())
    tools.sort()

    CMAP = "Purples_r"

    fig, axes = plt.subplots(2, 2, figsize=(5, 7))
    for col, bench in enumerate(
        ["Full", "strat-conf"]
    ):  # enumerate(df["Setting"].unique()):
        for row, refine in enumerate(df["Refined"].unique()):
            M = [[len(tools) for _ in truths] for _ in tools]
            M_annot = [["-" for _ in truths] for _ in tools]
            for hm_col, truth in enumerate(truths):
                df2 = df[
                    (df["Setting"] == bench)
                    & (df["ASM-Caller"] == truth)
                    & (df["Refined"] == refine)
                ]
                f1s = [(tool, f1) for tool, f1 in zip(df2["Caller"], df2["F1"])]
                f1s.sort(key=lambda x: x[1], reverse=True)
                # print(f1s)
                for rank, (tool, _) in enumerate(f1s, 1):
                    M[tools.index(tool)][hm_col] = rank
                    M_annot[tools.index(tool)][hm_col] = str(rank)
            g = sns.heatmap(
                M,
                ax=axes[row][col],
                annot=M_annot,
                fmt="",
                xticklabels=truths if row == 1 else False,
                yticklabels=tools if col == 0 else False,
                cbar=False,
                cmap=CMAP,
            )
            g.set_xticklabels(g.get_xticklabels(), rotation=30)
            if col == 0:
                g.set_ylabel("w/" + ("" if refine else "o") + " harmonization")
            else:
                g.set_ylabel("")
            if row == 0:
                title = "Full genome" if bench == "Full" else "Confident regions"
                axes[row][col].set_title(title)

    plt.suptitle(ref_name)
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "strip":
        main_strip()
    elif mode == "heat":
        main_heat()
    elif mode == "st":
        supp_table()
