import sys
from pysam import VariantFile


def main():
    vcf_fn = sys.argv[1]
    w = int(sys.argv[2])

    intervals = {}
    for record in VariantFile(vcf_fn):
        chrom, start, end = record.contig, record.start - 1000, record.stop + 1000

        gt = record.samples[0]["GT"]
        if chrom not in intervals:
            intervals[chrom] = [[start, end, [gt]]]
        else:
            if start > intervals[chrom][-1][1]:
                intervals[chrom].append([start, end, [gt]])
            else:
                intervals[chrom][-1] = [
                    min(intervals[chrom][-1][0], start),
                    max(intervals[chrom][-1][1], end),
                    intervals[chrom][-1][2] + [gt],
                ]
    for chrom, ints in intervals.items():
        for s, e, GTs in ints:
            if s < 0 or e < 0:
                continue
            first, second = False, False
            for gt in GTs:
                if gt[0] == 1:
                    first = True
                if gt[1] == 1:
                    second = True
            print(
                f"{chrom}:{s}-{e}",
                "first" if first else "",
                "second" if second else "",
                sep=" ",
            )


if __name__ == "__main__":
    main()
