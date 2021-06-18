
library(dplyr)
library(GenomicRanges)
library(purrr)
library(readr)

source("R/interval_simulation_functions.R")


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


qanbari_overlaps <- map(sweep_ranges,
                        function(x) subsetByOverlaps(x,
                                                     qanbari_ranges))

## Simulate overlaps

chroms <- read_tsv("data/galGal4.chrom.sizes",
                   col_names = FALSE)

autosomal_size <- sum(chroms$X2[chroms$X1 %in% paste("chr", 1:33, sep = "")])





sim1 <- simulate_overlaps(width(sweep_ranges[[1]]) - 1,
                          width(qanbari_ranges),
                          autosomal_size,
                          n_rep = 500)


sim2 <- simulate_overlaps(width(sweep_ranges[[2]]) - 1,
                          width(qanbari_ranges),
                          autosomal_size,
                          n_rep = 500)


sim3 <- simulate_overlaps(width(qanbari_ranges),
                          width(sweep_ranges[[3]]) - 1,
                          autosomal_size,
                          n_rep = 500)
    
print(mean(sim1))
print(quantile(sim1, 0.95))


print(mean(sim2))
print(quantile(sim2, 0.95))


print(mean(sim3))
print(quantile(sim3, 0.95))
