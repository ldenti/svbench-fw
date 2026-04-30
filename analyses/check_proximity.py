import sys
import os

from pysam import VariantFile
from intervaltree import IntervalTree

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme()

REFS = {"hg19": "GRCh37", "hg38": "GRCh38", "chm13": "T2T-CHM13"}
TRUTHS = {"dipcall": "dipcall", "svim-asm": "SVIM-asm", "hapdiff": "hapdiff"}
CALLERS = {
    "SVDSS-s4-q0": "SVDSS",
    "cutesv-s4-q20": "cuteSV",
    "debreak-s0-q20": "debreak",
    "sawfish-s0-q20": "sawfish",
    "severus-s4-q20": "severus",
    "sniffles-s0-q20": "sniffles",
    "svisionpro-s4-q20": "SVision-pro",
}


def get_uidx(record):
    return f"{record.contig}:{record.pos}:{record.ref}:{record.alts[0]}"


def vcf2trees(fn, chroms=[]):
    trees = {}
    for record in VariantFile(fn):
        if len(chroms) != 0 and record.contig not in chroms:
            continue
        chrom = record.contig
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        uidx = get_uidx(record)
        trees[chrom][record.pos : record.stop + 1] = uidx
    return trees


def vcf2trees_refine(fn):
    tp_trees, other_trees = {}, {}
    for record in VariantFile(fn):
        chrom = record.contig
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        if chrom not in tp_trees:
            tp_trees[chrom] = IntervalTree()
            other_trees[chrom] = IntervalTree()
        uidx = get_uidx(record)
        bd = record.samples[0]["BD"]

        if bd == "TP":
            tp_trees[chrom][record.pos : record.stop + 1] = uidx
        else:
            assert bd in ["FN", "FP"]
            other_trees[chrom][record.pos : record.stop + 1] = uidx
    return tp_trees, other_trees


def main():
    WD = sys.argv[1]

    # print(
    #     "Reference",
    #     "ASM-Caller",
    #     "Caller",
    #     "P_smp",
    #     "R_smp",
    #     "F1_smp",
    #     "P_cpx",
    #     "R_cpx",
    #     "F1_cpx",
    #     sep=",",
    # )

    df = []
    for ref in REFS:
        for truth in TRUTHS:

            fp = os.path.join(
                WD, ref, "asmcallsets-giab", f"{truth}.noinfo.regions-w500.txt"
            )
            regions = []
            for line in open(fp):
                chrom, se = line.split(" ")[0].split(":")
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                s, e = se.split("-")
                s, e = int(s), int(e)
                regions.append((chrom, s, e))

            for caller in CALLERS:
                wd = os.path.join(WD, ref, "truvari", "giab", truth, caller)

                print(ref, truth, caller, file=sys.stderr)

                # if False:
                #     # no refine
                #     basetp_trees = vcf2trees(wd + "/tp-base.vcf.gz")
                #     comptp_trees = vcf2trees(wd + "/tp-comp.vcf.gz")
                #     fp_trees = vcf2trees(wd + "/fp.vcf.gz")
                #     fn_trees = vcf2trees(wd + "/fn.vcf.gz")
                # refine
                basetp_trees, fn_trees = vcf2trees_refine(
                    wd + "/ga4gh_with_refine.base.vcf.gz"
                )
                comptp_trees, fp_trees = vcf2trees_refine(
                    wd + "/ga4gh_with_refine.comp.vcf.gz"
                )

                smp = [0, 0, 0, 0]
                cpx = [0, 0, 0, 0]

                for chrom, start, end in regions:
                    basetps, fns, comptps, fps = set(), set(), set(), set()
                    if chrom in basetp_trees:
                        basetps = basetp_trees[chrom].overlap(start, end)
                    if chrom in fn_trees:
                        fns = fn_trees[chrom].overlap(start, end)
                    if chrom in comptp_trees:
                        comptps = comptp_trees[chrom].overlap(start, end)
                    if chrom in fp_trees:
                        fps = fp_trees[chrom].overlap(start, end)

                    # if basetps + fns == 0:
                    #     # XXX: do we want this for precision?
                    #     continue

                    if len(basetps) + len(fns) + len(fps) == 0:
                        continue

                    if len(basetps) + len(fns) == 1:
                        smp[0] += len(basetps)
                        smp[1] += len(fns)
                        smp[2] += len(comptps)
                        smp[3] += len(fps)
                    elif len(basetps) + len(fns) > 1:
                        cpx[0] += len(basetps)
                        cpx[1] += len(fns)
                        cpx[2] += len(comptps)
                        cpx[3] += len(fps)

                    # region = "Simple" if basetps + fns <= 1 else "Multi"
                    # df.append(
                    #     [
                    #         f"{c}:{start}-{end}",
                    #         region,
                    #         basetps,
                    #         fns,
                    #         comptps,
                    #         fps,
                    #     ]
                    # )

                print(smp, cpx, file=sys.stderr)
                p_smp = smp[2] / (smp[2] + smp[3])
                r_smp = smp[0] / (smp[0] + smp[1])
                f_smp = 2 * p_smp * r_smp / (p_smp + r_smp)

                p_cpx = cpx[2] / (cpx[2] + cpx[3])
                r_cpx = cpx[0] / (cpx[0] + cpx[1])
                f_cpx = 2 * p_cpx * r_cpx / (p_cpx + r_cpx)

                p = (smp[2] + cpx[2]) / (smp[2] + cpx[2] + smp[3] + cpx[3])
                r = (smp[0] + cpx[0]) / (smp[0] + cpx[0] + smp[1] + cpx[1])
                f = 2 * p * r / (p + r)

                # print(
                #     ref,
                #     truth,
                #     caller,
                #     p_smp,
                #     r_smp,
                #     f_smp,
                #     p_cpx,
                #     r_cpx,
                #     f_cpx,
                #     p,
                #     r,
                #     f,
                #     sep=",",
                # )

                # df.append([ref, truth, caller, f_smp, f_cpx])
                df.append([ref, truth, caller, "Singleton", p_smp, r_smp, f_smp] + smp)
                df.append([ref, truth, caller, "Multi", p_cpx, r_cpx, f_cpx] + cpx)

    # df = pd.DataFrame(
    #     df, columns=["Reference", "Truth", "Caller", "F1 (simple)", "F1 (complex)"]
    # )

    df = pd.DataFrame(
        df,
        columns=[
            "Reference",
            "Truth",
            "Caller",
            "Region",
            "P",
            "R",
            "F1",
            "TP-base",
            "FN",
            "TP-comp",
            "FP",
        ],
    )
    print(df.to_csv(index=False))


def plot_f1():
    csv_fn = sys.argv[1]
    df = pd.read_csv(csv_fn)

    df["Reference"] = df["Reference"].replace(REFS)
    df["Truth"] = df["Truth"].replace(TRUTHS)
    df["Caller"] = df["Caller"].replace(CALLERS)

    fig, axes = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(9, 8))
    for col, ref in enumerate(REFS):
        ref = REFS[ref]
        for row, truth in enumerate(TRUTHS):
            truth = TRUTHS[truth]
            ax = axes[row][col]
            subdf = df[(df["Reference"] == ref) & (df["Truth"] == truth)]
            sns.barplot(
                data=subdf,
                x="Caller",
                y="F1",
                hue="Region",
                legend=True if col == 0 and row == 0 else False,
                ax=ax,
            )

            if row == 0:
                ax.set_title(ref)
            if row == 2:
                ax.set_xlabel("Precision")
                ax.set_xticklabels(ax.get_xticklabels(), rotation=60)

            if col == 0:
                ax.set_ylabel(f"{truth}\nRecall")

    legend = axes[0][0].get_legend()  # Legend instance
    legend.remove()
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles=handles,
        labels=labels,
        loc="center right",
        frameon=False,
    )

    plt.tight_layout()
    fig.subplots_adjust(right=0.83)  # make room for the legend
    plt.show()


def plot_pr():
    csv_fn = sys.argv[1]
    pdf_fn = sys.argv[2]

    df = pd.read_csv(csv_fn)

    df["Reference"] = df["Reference"].replace(REFS)
    df["Truth"] = df["Truth"].replace(TRUTHS)
    df["Caller"] = df["Caller"].replace(CALLERS)

    fig, axes = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(9, 8))
    for col, ref in enumerate(REFS):
        ref = REFS[ref]
        for row, truth in enumerate(TRUTHS):
            truth = TRUTHS[truth]
            ax = axes[row][col]
            subdf = df[(df["Reference"] == ref) & (df["Truth"] == truth)]
            sns.scatterplot(
                data=subdf,
                x="P",
                y="R",
                hue="Caller",
                style="Region",
                ax=ax,
                legend=True if col == 0 and row == 0 else False,
            )
            ax.set_xlim(0.1, 1.0)
            ax.set_ylim(0.1, 1.0)

            if row == 0:
                ax.set_title(ref)
            if row == 2:
                ax.set_xlabel("Precision")
            if col == 0:
                ax.set_ylabel(f"{truth}\nRecall")

            # F1 iso-curves
            f_scores = np.linspace(0.2, 0.9, 8)  # F1 levels to draw
            recall_grid = np.linspace(0.01, 1.0, 500)
            precision_grid = np.linspace(0.01, 1.0, 500)
            R, P = np.meshgrid(recall_grid, precision_grid)
            F1 = 2 * (P * R) / (P + R)
            cs = ax.contour(
                R,
                P,
                F1,
                levels=f_scores,
                colors="0.6",
                linestyles="dashed",
                linewidths=0.8,
            )
            plt.clabel(cs, inline=1, fmt="F1=%.2f", fontsize=8)

            ax.set_xticks([0.25, 0.5, 0.75, 1.0])
            ax.set_yticks([0.25, 0.5, 0.75, 1.0])

    legend = axes[0][0].get_legend()  # Legend instance
    legend.remove()
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles=handles,
        labels=labels,
        loc="center right",
        frameon=False,
    )

    plt.tight_layout()
    fig.subplots_adjust(right=0.83)  # make room for the legend
    # plt.show()
    plt.savefig(pdf_fn)


def plot2():
    csv_fn = sys.argv[1]
    df = pd.read_csv(csv_fn)

    df["Reference"] = df["Reference"].replace(REFS)
    df["Truth"] = df["Truth"].replace(TRUTHS)
    df["Caller"] = df["Caller"].replace(CALLERS)

    df["FN"] = df["TP-base"] + df["FN"]
    df["FP"] = df["TP-comp"] + df["FP"]

    df = df.melt(
        id_vars=["Reference", "Truth", "Caller", "Region"],
        # value_vars=["TP-base", "FN", "Total-base", "TP-comp", "FP", "Total-comp"],
        value_vars=["TP-base", "FN", "FP", "TP-comp"],
        var_name="Class",
        value_name="Count",
    )
    print(df)

    fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(9, 8))
    ref = "T2T-CHM13"
    for row, truth in enumerate(TRUTHS):
        truth = TRUTHS[truth]
        subdf = df[
            (df["Reference"] == ref)
            & (df["Truth"] == truth)
            & (df["Region"] == "Singleton")
        ]

        sns.barplot(
            data=subdf[df["Class"].isin(["FN", "FP"])],
            # data=subdf[df["Class"].isin(["Total-base", "Total-comp"])],
            x="Caller",
            y="Count",
            hue="Class",
            alpha=0.5,
            legend=True if row == 0 else False,
            ax=axes[row][0],
        )
        sns.barplot(
            data=subdf[df["Class"].isin(["TP-base", "TP-comp"])],
            x="Caller",
            y="Count",
            hue="Class",
            legend=True if row == 0 else False,
            ax=axes[row][0],
        )

        subdf = df[
            (df["Reference"] == ref)
            & (df["Truth"] == truth)
            & (df["Region"] == "Multi")
        ]
        sns.barplot(
            data=subdf[df["Class"].isin(["FN", "FP"])],
            # data=subdf[df["Class"].isin(["Total-base", "Total-comp"])],
            x="Caller",
            y="Count",
            hue="Class",
            alpha=0.5,
            legend=False,
            ax=axes[row][1],
        )
        sns.barplot(
            data=subdf[df["Class"].isin(["TP-base", "TP-comp"])],
            x="Caller",
            y="Count",
            hue="Class",
            legend=False,
            ax=axes[row][1],
        )

        # sns.barplot(data=subdf, x="Caller", y="Count", hue="Class", ax=axes[row][1])

        if row == 0:
            axes[0][0].set_title("Singleton")
            axes[0][1].set_title("Multi")
        if row == 2:
            axes[row][0].set_xticklabels(axes[row][0].get_xticklabels(), rotation=30)
            axes[row][1].set_xticklabels(axes[row][0].get_xticklabels(), rotation=30)
        axes[row][0].set_ylabel(truth)

    legend = axes[0][0].get_legend()  # Legend instance
    legend.remove()
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles=handles,
        labels=labels,
        ncol=4,
        loc="lower center",
        frameon=False,
    )

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.17)  # make room for the legend
    plt.show()


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "pr":
        plot_pr()
    elif mode == "f1":
        plot_f1()
    elif mode == "plot2":
        plot2()
    else:
        main()
