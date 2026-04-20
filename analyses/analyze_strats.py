import argparse
import os
import numpy as np

STRAT_MAP = {
    "easy": "Easy",
    "lowcomplex": "Low Complexity.",
    "lowmap": "Low Mappability",
    "segdups": "Segmental Duplications",
}

GENOME_L = (
    3117292070  # from $(cut -f2 /projects6/svdss2/chm13v2.0.fa.fai | paste -sd+ | bc)
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("WD")

    args = parser.parse_args()

    print(
        "Reference",
        "Region",
        "Count",
        "Min",
        "Q25",
        "Q50",
        "Q75",
        "Max",
        "Avg",
        "Std",
        "Sum",
        "Coverage",
        sep=",",
    )
    for reference in ["chm13", "hg38", "hg19"]:
        for strat in ["easy", "lowcomplex", "lowmap", "segdups"]:
            data = []
            bed_fn = os.path.join(args.WD, "input", "strats", strat, f"{reference}.bed")
            for line in open(bed_fn):
                chrom, s, e, *_ = line.strip("\n").split("\t")
                s, e = int(s), int(e)
                l = e - s
                data.append(l)
            strat = STRAT_MAP[strat]
            data = np.array(data)
            print(
                reference,
                strat,
                len(data),
                data.min(),
                np.quantile(data, 0.25),
                np.median(data),
                np.quantile(data, 0.75),
                data.max(),
                round(np.mean(data), 3),
                round(np.std(data), 3),
                sum(data),
                round(sum(data) / GENOME_L, 3),
                sep=",",
            )


if __name__ == "__main__":
    main()
