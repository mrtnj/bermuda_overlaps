
library(dplyr)
library(GenomicRanges)
library(purrr)
library(readr)



## Bermuda sweeps

files <- system("ls data/Bermuda*_regions",
                intern = TRUE)

sweeps <- map(files, read_tsv)

sweep_ranges <- map(sweeps, makeGRangesFromDataFrame)

names(sweep_ranges) <- c("Hp", "Tajima's D", "Both Hp and Tajima's D")



## Qanbari sweeps

qanbari <- read_tsv("data/pgen.1007989.s009_Qanbari_table.txt")

qanbari$id <- paste(qanbari$CHR,
                    qanbari$BIN_START,
                    qanbari$BIN_END,
                    qanbari$Contrast,
                    sep = "_")


## Prepare for lifting

qanbari_bed <- data.frame(chr = paste("chr", qanbari$CHR, sep = ""),
                          start = qanbari$BIN_START - 1,
                          end = qanbari$BIN_END,
                          name = qanbari$id)

options(scipen = 1e9)

write.table(qanbari_bed,
            file = "data/qanbari_to_lift.bed",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


## Lifted to Galgal4

qanbari_lift <- read_tsv("data/qanbari_hglft_galgal4_0.7.bed",
                         col_names = FALSE)
                         


qanbari_lifted <- inner_join(qanbari_lift, qanbari, by = c("X4" = "id"))


qanbari_ranges <- GRanges(seqnames = qanbari_lifted$X1,
                          ranges = IRanges(qanbari_lifted$X2 + 1,
                                           qanbari_lifted$X3),
                          mcols = qanbari_lifted)


findOverlaps(qanbari_ranges, sweep_ranges[[1]])
