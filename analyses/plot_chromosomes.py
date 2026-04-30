import sys
import glob
from pysam import VariantFile
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


def get_uidx(v):
    return f"{v.contig}:{v.pos}:{v.ref}:{v.alts[0]}"


def main():
    vcf_dir = sys.argv[1]
    region = sys.argv[2]  # region (should work for single chromosome too)

    chrom = region
    start, end = -1, -1
    if ":" in region:
        chrom = region.split(":")[0]
        start, end = region.split(":")[1].split("-")
        start = int(start)
        end = int(end)

    if start != -1:
        chrom_l = end - start + 1
    else:
        vcf_fn = f"{vcf_dir}/dipcall.vcf.gz"  # assuming dipcall to be there
        vcf = VariantFile(vcf_fn)
        for name in vcf.header.contigs:
            if name != chrom:
                continue
            chrom_l = vcf.header.contigs[name].length
            break
        vcf.close()
    print(chrom_l)
    assert chrom_l != 0

    start = 0 if start == -1 else start
    end = chrom_l if end == -1 else end

    TPs = {"dipcall": [], "svim-asm": [], "hapdiff": []}
    truths = list(TPs.keys())
    print(truths)

    data = []

    for truth in truths:
        vcf_fn = f"{vcf_dir}/{truth}.vcf.gz"
        for record in VariantFile(vcf_fn).fetch(region=region):
            l = len(record.alts[0]) - len(record.ref)
            if abs(l) < 50:
                continue
            # if record.pos < start or record.pos > end:
            #     continue
            gt1, gt2 = record.samples[0]["GT"]
            gt1 = gt1 if gt1 != None else 0
            gt2 = gt2 if gt2 != None else 0
            svtype = "INS" if l > 0 else "DEL"
            l = 1 if svtype == "INS" else l

            color = "green"
            if gt1 == 0 and gt2 == 0:
                continue
            elif gt1 > 0 and gt2 > 0:
                color = "red"
            if gt1 > 0:
                data.append(
                    [
                        truth,
                        record.contig,
                        record.pos,
                        record.pos + l,
                        l,
                        svtype,
                        color,
                        1,
                    ]
                )
            if gt2 > 0:
                data.append(
                    [
                        truth,
                        record.contig,
                        record.pos,
                        record.pos + l,
                        l,
                        svtype,
                        color,
                        2,
                    ]
                )

    regions = pd.DataFrame(
        data,
        columns=["truth", "chrom", "start", "end", "width", "type", "color", "haplo"],
    )
    print(regions)
    for truth in TPs:
        print(f"=== {truth} ===")
        print(len(regions[(regions["truth"] == truth)]))
        print(len(regions[(regions["truth"] == truth) & (regions["haplo"] == 1)]))
        print(len(regions[(regions["truth"] == truth) & (regions["haplo"] == 2)]))
        for color in ["green", "red"]:
            print(
                color,
                len(regions[(regions["truth"] == truth) & (regions["color"] == color)]),
            )

    figsize = (8, 6)
    # sep=10

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # add regions to the axis

    chrom_height = 1  # Height of each ideogram
    chrom_spacing = 1  # Spacing between consecutive ideograms

    gene_height = 0.4  # Height of the gene track < `chrom_spacing`
    gene_padding = 0.1  # Padding between the top of a gene track ideogram

    # Decide which chromosomes to use
    # chromosome_list = list(chromosomes)

    # Keep track of the y positions genes for each chromosome
    # and center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}
    for truth in TPs:
        tt = truth + ".1"
        chrom_ybase[tt] = ybase
        chrom_centers[truth] = ybase + chrom_height + chrom_spacing / 8
        gene_ybase[tt] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing / 4

        tt = truth + ".2"
        chrom_ybase[tt] = ybase
        gene_ybase[tt] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # add chromsizes to ax
    ideo = []
    for truth in TPs:
        ideo.append([truth + ".1", end, "whitesmoke", start, end - start + 1])
        ideo.append([truth + ".2", end, "whitesmoke", start, end - start + 1])
    ideo = pd.DataFrame(ideo, columns=["truth", "end", "colors", "start", "width"])

    for truth, group in ideo.groupby("truth"):
        yrange = (chrom_ybase[truth], chrom_height)
        xranges = group[["start", "width"]].values
        ax.broken_barh(xranges, yrange, facecolors=group["colors"], edgecolor="black")

    # # legend
    # #leg = []
    # #leg_lab = []

    for i, r in enumerate(regions):
        # add legend element per regions
        # leg_lab.append(f"regions{i+1}")
        # leg.append(Line2D([0], [0], color=color, lw=4))
        # chromosome_collections(ax, regions, chrom_ybase, chrom_height, edgecolor=color)

        for h in [1, 2]:
            for truth, group in regions[regions["haplo"] == h].groupby("truth"):
                yrange = (chrom_ybase[truth + f".{h}"], chrom_height)
                xranges = group[["start", "width"]].values
                ax.broken_barh(
                    xranges,
                    yrange,
                    facecolors=group["color"],
                    edgecolor=group["color"],
                    alpha=0.1,
                )

    # add to ax
    ax.set_yticks([chrom_centers[i] for i in truths])
    ax.set_yticklabels(truths)

    # ax.legend(leg, leg_lab, loc=4)
    ax.axis("tight")

    # # ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))

    # labels
    ax.set_title(f"{chrom}", fontsize=14)
    plt.xlabel("Chromosome Position (Mbp)", fontsize=14)
    plt.ylabel("Truth", fontsize=14)

    plt.tight_layout()
    # if savefig is True:
    #     plt.savefig(f'{species}_karyopype.png')
    plt.show()


if __name__ == "__main__":
    main()
