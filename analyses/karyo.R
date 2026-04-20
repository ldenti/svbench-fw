library(karyoploteR)
## library(GenomicRanges)

## library(BSgenome)
## available.genomes()

args <- commandArgs(trailingOnly = TRUE)

fn1 <- args[1]
fn2 <- args[2]
ref <- args[3]
title <- args[4]
pdf_fn <- args[5]

## fn1 <- "/data/svbench-natcomm/HG002/chm13/asmcallsets-hprc/dipcall.merged.bed"
## fn2 <- "/data/svbench-natcomm/HG002/chm13/asmcallsets-hprc/hapdiff.merged.bed"
## ref <- "chm13"

## fn1 <- "/data/svbench-natcomm/hg002/giab/dipcall.merged.bed"
dipcall <- read.table(fn1, header = FALSE)
colnames(dipcall) <- c("chrom", "start", "end")
dipcall$chrom <- gsub("^chr", "", dipcall$chrom)

## fn2 <- "/data/svbench-natcomm/hg002/giab/confident_regions.merged.bed"
hapdiff <- read.table(fn2, header = FALSE)
colnames(hapdiff) <- c("chrom", "start", "end")
hapdiff$chrom <- gsub("^chr", "", hapdiff$chrom)

## # Define the chromosomes you want to visualize
## chromosomes <- c("chr21")

## ## # Initialize the Karyoploter
## ## kp <- karyoploteR::initKaryotype()

## ## # Add data to the plot
## ## kpPlotRegions(kp, data = bed_data, col = "red")

## kp <- plotKaryotype(genome="hg19", chromosomes=c("chrX", "chrY"))


# 1. Create the base plot
pp <- getDefaultPlotParams(plot.type = 1)
pp$data1height <- 50

genome <- ""
if (ref == "chm13") {
  genome <- "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
} else if (ref == "hg38") {
  genome <- "BSgenome.Hsapiens.NCBI.GRCh38"
} else if (ref == "hg19") {
  genome <- "BSgenome.Hsapiens.1000genomes.hs37d5"
}

print(genome)
print(ref)

chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
pdf(pdf_fn, width = 7, height = 7)

kp <- plotKaryotype(genome = genome, chromosomes = chroms, plot.type = 1, plot.params = pp)
# , chromosomes="21")

# 2. Define vertical positions for 3 panels above (r0=bottom, r1=top)
# The ideogram is usually near y=0, we move up.
panel1 <- autotrack(current.track = 1, total.tracks = 2, margin = 0.1)
panel2 <- autotrack(current.track = 2, total.tracks = 2, margin = 0.1)

# panel3 <- autotrack(current.track=3, total.tracks=3, margin=0.1)

# 3. Add data to the 3 panels

kpPlotRegions(kp,
  data = dipcall,
  data.panel = 1, r0 = 0, r1 = 0.5, # r0 = panel1$r0, r1 = panel1$r1,
  col = "#DDDD00", border = "#000000"
)

kpPlotRegions(kp,
  data = hapdiff,
  data.panel = 2, r0 = 0.5, r1 = 1, # r0 = panel2$r0, r1 = panel2$r1,
  col = "#00DDDD", border = "#000000"
)

kpAddMainTitle(kp, title)

legend("bottomright", legend = c("dipcall", "hapdff"), fill = c("#DDDD00", "#00DDDD"), bty = "n")

dev.off()

## # Panel 1 (closest to ideogram)
## kpPoints(kp, chr="chr21", x=10e6, y=0.5, r0=panel1$r0, r1=panel1$r1, col="red")
## # Panel 2
## kpLines(kp, chr="chr21", x=c(20e6, 30e6), y=c(0.2, 0.8), r0=panel2$r0, r1=panel2$r1, col="blue")
## # Panel 3 (furthest from ideogram)
## kpRect(kp, chr="chr21", x0=40e6, x1=50e6, y0=0, y1=1, r0=panel3$r0, r1=panel3$r1, col="green")
