import sys
import glob
import argparse
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=2)
sns.set(style="whitegrid")

REFSEQS = ["GRCh37", "GRCh38", "T2T-CHM13"]
TRUTHS = ["dipcall", "svim-asm", "hapdiff"]


def parse_dir(ddir, asm, confident=False, refseq=""):
    # we may not need this dict and do everything on df but it's more convenient to me
    truths = {}
    data = []
    conf = ".confident" if confident else ""
    for vcf_fn in glob.glob(f"{ddir}/asmcallsets-{asm}{conf}/*.vcf.gz"):
        if "noinfo" in vcf_fn:
            continue
        name = vcf_fn.split("/")[-1].split(".")[0]
        if name not in TRUTHS:
            continue
        truths[name] = {}
        for record in VariantFile(vcf_fn):
            l = len(record.alts[0]) - len(record.ref)
            # filters=[x.name for x in record.filter.values()]
            # if name == "dipcall" and len(filters) > 0 and ("GAP1" in filters or "GAP2" in filters):
            #     continue
            # elif name != "dipcall" and "PASS" not in filters:
            #     continue

            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0

            if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
                data.append(
                    [
                        refseq,
                        name,
                        record.contig,
                        record.pos,
                        l,
                        "INS" if l > 0 else "DEL",
                        gt1 == gt2,
                    ]
                )
                if record.contig not in truths[name]:
                    truths[name][record.contig] = []
                truths[name][record.contig].append(record.pos)
    return data, truths


def parse_pafs(ddir, asm, confident=False, refseq=""):
    data = []
    conf = ".confident" if confident else ""
    for paf_fn in glob.glob(f"{ddir}/asmcallsets-{asm}{conf}/*.haps-w*.paf"):
        name = paf_fn.split("/")[-1].split(".")[0]
        if name not in TRUTHS:
            continue
        for line in open(paf_fn):
            line = line.strip("\n").split("\t")
            nm = int(line[12].split(":")[-1])
            tp = line[16][-1]
            de = float(line[20].split(":")[-1])
            if tp == "P":
                data.append([refseq, name, nm, de])
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--confident", action="store_true")
    parser.add_argument("--asm", required=True, type=str)
    parser.add_argument("--dist", default=500)
    parser.add_argument("-o", default="")
    parser.add_argument("WD")

    args = parser.parse_args()

    df = []
    d, t2t_truth = parse_dir(args.WD + "/chm13", args.asm, args.confident, "T2T-CHM13")
    df += d
    d, hg38_truth = parse_dir(args.WD + "/hg38", args.asm, args.confident, "GRCh38")
    df += d
    d, hg19_truth = parse_dir(args.WD + "/hg19", args.asm, args.confident, "GRCh37")
    df += d

    df = pd.DataFrame(
        df, columns=["RefSeq", "Truth", "Chrom", "Pos", "Len", "Type", "GT"]
    )

    fig, axes = plt.subplots(4, 3, sharey="row", figsize=(9, 11))

    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        df2 = subdf.groupby(["Truth", "Type"]).count()
        sns.barplot(
            data=df2,
            x="Truth",
            order=TRUTHS,
            y="Chrom",
            hue="Type",
            legend=True if i == 2 else None,
            palette="Set2",
            ax=axes[0][i],
        )
        axes[0][i].set_xlabel("")  # Truth
        axes[0][i].tick_params(axis="x", labelrotation=0)
        # axes[0][i].set_ylim(0, 30000)  # Count
        axes[0][i].set_ylabel("")  # Count
        if i == 0:
            axes[0][i].set_ylabel("(a)\nCount")
        if i == 2:
            # move legends
            axes[0][i].legend(ncols=2, handletextpad=0.2, columnspacing=0.2)
            # sns.move_legend(axes[0][i], "center left", bbox_to_anchor=(1, 0.5))

        # ax1.bar_label(ax1.containers[0])
        # ax1.bar_label(ax1.containers[1])
        axes[0][i].set_title(refseq)

        # variant length distribution per truthset
        sns.histplot(
            data=subdf[abs(subdf["Len"]) <= 500],
            x="Len",
            hue="Truth",
            hue_order=TRUTHS,
            element="poly",
            fill=False,
            legend=False,  # True if i == 2 else None,
            ax=axes[1][i],
        )
        axes[1][i].set_xlabel("Length")
        axes[1][i].set_xticks([-500, -250, 0, 250, 500])
        axes[1][i].set_ylabel("")  # Count
        # axes[1][i].set_ylim(0, 2500)
        if i == 0:
            axes[1][i].set_ylabel("(b)\nCount")
        # if i == 2:
        #     sns.move_legend(axes[1][i], "center left", bbox_to_anchor=(1, 0.5))

    # neighbor distribution per truthset
    print("=== Statistics on neighbouring calls")
    for i, truths in enumerate([hg19_truth, hg38_truth, t2t_truth]):
        df2 = []
        for truth in truths:
            neighbors = []
            for chrom in truths[truth]:
                last_p = truths[truth][chrom][0]
                neighbors.append(0)
                for p in truths[truth][chrom][1:]:
                    if p - last_p > args.dist:
                        neighbors.append(0)
                    else:
                        neighbors[-1] += 1
                    last_p = p
            d = {}
            for x in neighbors:
                k = str(x)
                if x >= 2:
                    k = "2+"
                d[k] = d[k] + 1 if k in d else 1
            for k, v in d.items():
                df2.append([truth, k, v])
            refseq = "T2T"
            if i == 0:
                refseq = "hg19"
            elif i == 1:
                refseq = "hg38"
            print(
                refseq,
                truth,
                d["0"],
                d["1"],
                d["2+"],
                (d["1"] + d["2+"]) / (d["0"] + d["1"] + d["2+"]),
            )

        df2 = pd.DataFrame(df2, columns=["Truth", f"#Neighbors-{args.dist}bp", "Count"])

        sns.barplot(
            data=df2,
            x=f"#Neighbors-{args.dist}bp",
            order=["0", "1", "2+"],
            y="Count",
            hue="Truth",
            hue_order=TRUTHS,
            ax=axes[2][i],
            legend=False,  # True if i == 2 else None,
        )
        # axes[2][i].set_ylim(0, 30000)
        axes[2][i].set_ylabel("")  # Count
        if i == 0:
            axes[2][i].set_ylabel("(c)\nCount")
        # if i == 2:
        #     # move legends
        #     sns.move_legend(axes[2][i], "center left", bbox_to_anchor=(1, 0.5))
    print("====================================")

    # --- NM ttmars-like
    df = []
    df += parse_pafs(args.WD + "/chm13", args.asm, args.confident, "T2T-CHM13")
    df += parse_pafs(args.WD + "/hg38", args.asm, args.confident, "GRCh38")
    df += parse_pafs(args.WD + "/hg19", args.asm, args.confident, "GRCh37")

    df = pd.DataFrame(df, columns=["RefSeq", "Truth", "NM", "de"])

    for i, refseq in enumerate(REFSEQS):
        subdf = df[df["RefSeq"] == refseq]
        for truth in TRUTHS:
            print("NM", refseq, truth)
            print(subdf[subdf["Truth"] == truth]["NM"].describe())
        # sns.histplot(
        #     data=subdf,
        #     x="NM",
        #     hue="Truth",
        #     hue_order=TRUTHS,
        #     binrange=[0, 50],
        #     discrete=True,
        #     element="step",
        #     legend=True if i == 2 else None,
        #     ax=axes[3][i],
        # )
        sns.boxplot(
            data=subdf,  # [subdf["NM"] < 100],
            x="Truth",
            y="NM",
            hue="Truth",
            order=TRUTHS,
            hue_order=TRUTHS,
            showfliers=False,
            ax=axes[3][i],
        )
        # axes[3][i].set_ylim(-10, 325)
        if i == 0:
            axes[3][i].set_ylabel("(d)\nNM")
        else:
            axes[3][i].set_ylabel("")
        # if i == 2:
        #     # move legends
        #     sns.move_legend(axes[3][i], "center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
