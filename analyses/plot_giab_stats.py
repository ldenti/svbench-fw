import sys
import glob
import argparse
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=2)
sns.set(style="whitegrid")


def parse_vcf(vcf_fn, tag):
    data = []
    for record in VariantFile(vcf_fn):
        l = len(record.alts[0]) - len(record.ref)
        # filters=[x.name for x in record.filter.values()]
        # if name == "dipcall" and len(filters) > 0 and ("GAP1" in filters or "GAP2" in filters):
        #     continue
        # elif name != "dipcall" and "PASS" not in filters:
        #     continue

        gt1, gt2 = 0, 0
        if len(record.samples[0]["GT"]) == 2:
            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0
        else:
            gt1 = record.samples[0]["GT"][0]
            gt1 = gt1 if gt1 != None else 0

        if abs(l) >= 50 and (gt1 > 0 or gt2 > 0):
            data.append(
                [
                    tag,
                    l,
                    "INS" if l > 0 else "DEL",
                ]
            )
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", default="")
    parser.add_argument("WD")

    args = parser.parse_args()

    df = []
    df += parse_vcf(args.WD + "/input/giab06/hg19.vcf.gz", "GRCh37 (v0.6)")
    df += parse_vcf(args.WD + "/input/giab50/hg19.vcf.gz", "GRCh37 (v5.0q)")
    df += parse_vcf(args.WD + "/input/giab50/hg38.vcf.gz", "GRCh38 (v5.0q)")
    df += parse_vcf(args.WD + "/input/giab50/chm13.vcf.gz", "CHM13 (v5.0q)")

    df_conf = []
    df_conf += parse_vcf(args.WD + "/input/giab06/hg19-conf.vcf.gz", "GRCh37 (v0.6)")
    df_conf += parse_vcf(args.WD + "/input/giab50/hg19-conf.vcf.gz", "GRCh37 (v5.0q)")
    df_conf += parse_vcf(args.WD + "/input/giab50/hg38-conf.vcf.gz", "GRCh38 (v5.0q)")
    df_conf += parse_vcf(args.WD + "/input/giab50/chm13-conf.vcf.gz", "CHM13 (v5.0q)")

    ORDER = ["GRCh37 (v0.6)", "GRCh37 (v5.0q)", "GRCh38 (v5.0q)", "CHM13 (v5.0q)"]

    df = pd.DataFrame(df, columns=["Callset", "Len", "Type"])
    df_conf = pd.DataFrame(df_conf, columns=["Callset", "Len", "Type"])

    fig, axes = plt.subplots(2, 2, sharey="row", figsize=(11, 11))

    sns.barplot(
        data=df.groupby(["Callset", "Type"]).count(),
        x="Callset",
        order=ORDER,
        y="Len",
        hue="Type",
        legend=False,
        palette="Set2",
        ax=axes[0][0],
    )
    sns.barplot(
        data=df_conf.groupby(["Callset", "Type"]).count(),
        x="Callset",
        order=ORDER,
        y="Len",
        hue="Type",
        legend=True,
        palette="Set2",
        ax=axes[0][1],
    )

    axes[0][0].set_title("Full genome")
    axes[0][1].set_title("Confident regions")

    axes[0][0].set_xlabel("")
    axes[0][0].set_ylabel("(a)\nCount")
    axes[0][1].set_xlabel("")
    axes[0][0].set_xticklabels(axes[0][0].get_xticklabels(), rotation=10)
    axes[0][1].set_xticklabels(axes[0][1].get_xticklabels(), rotation=10)

    sns.histplot(
        data=df[abs(df["Len"]) <= 500],
        x="Len",
        hue="Callset",
        hue_order=ORDER,
        element="poly",
        fill=False,
        legend=False,
        ax=axes[1][0],
    )
    sns.histplot(
        data=df_conf[abs(df_conf["Len"]) <= 500],
        x="Len",
        hue="Callset",
        hue_order=ORDER,
        element="poly",
        fill=False,
        legend=True,
        ax=axes[1][1],
    )
    axes[1][0].set_ylabel("(b)\nCount")

    plt.tight_layout()
    if args.o == "":
        plt.show()
    else:
        plt.savefig(args.o)


if __name__ == "__main__":
    main()
