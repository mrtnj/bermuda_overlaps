
library(dplyr)
library(GenomicRanges)
library(purrr)
library(readr)

source("R/interval_simulation_functions.R")


## Bermuda sweeps

files <- c("data/ehh_sweeps.csv", "data/clr_sweeps.csv", "data/tajima_sweeps.csv")

sweeps <- map(files, read_csv)

names(sweeps) <- c("EHH", "CLR", "Tajima")

sweeps_split <- map(sweeps, function(x) split(x, x$population))

sweeps_list <- list_flatten(sweeps_split)

sweep_ranges <- map(sweeps_list, makeGRangesFromDataFrame)


consensus_files <- c("data/consensus_bermuda.csv", "data/consensus_kauai.csv")

consensus <- map(consensus_files, read_csv)

names(consensus) <- c("bermuda", "kauai")

consensus_ranges <- map(
    consensus,
    function(x) {
        GenomicRanges::reduce(GRanges(seqnames = x$chr_name,
                                      ranges = IRanges(x$region_start, x$region_end)))
    }
)



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


## Lifted to Galgal6

qanbari_lift <- read_tsv("data/qnbari_hglft_galgal6.bed",
                         col_names = FALSE)
                         


qanbari_lifted <- inner_join(qanbari_lift, qanbari, by = c("X4" = "id"))

qanbari_lifted$X1 <- sub(qanbari_lifted$X1,
                         pattern = "chr",
                         replacement = "")


qanbari_ranges <- GenomicRanges::reduce(GRanges(seqnames = qanbari_lifted$X1,
                                                ranges = IRanges(qanbari_lifted$X2 + 1,
                                                                 qanbari_lifted$X3),
                                                mcols = qanbari_lifted))


qanbari_overlaps <- map(sweep_ranges,
                        function(x) subsetByOverlaps(x,
                                                     qanbari_ranges))

## Simulate overlaps

chroms <- read_tsv("data/galGal6.chrom.sizes",
                   col_names = FALSE)

autosomal_size <- sum(chroms$X2[chroms$X1 %in% paste("chr", 1:33, sep = "")])



analyse_sim <- function(sim, n) {
    
    q <- quantile(sim, 0.95)
    pt <- (sum(sim >= n) + 1)/1000
    
    tibble(q95 = q,
           p_less_than = pt)
    
}



## EHH

sim_bermuda_kauai_ehh <- simulate_overlaps(width(sweep_ranges$EHH_bermuda) - 1,
                                           width(sweep_ranges$EHH_kauaii) - 1,
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_bermuda_kauai_ehh, 17))

sim_qanbari_bermuda_ehh <- simulate_overlaps(width(sweep_ranges$EHH_bermuda) - 1,
                                             width(qanbari_ranges),
                                             autosomal_size,
                                             n_rep = 1000)

print(analyse_sim(sim_qanbari_bermuda_ehh, 5))

sim_qanbari_kauai_ehh <- simulate_overlaps(width(sweep_ranges$EHH_kauai) - 1,
                                           width(qanbari_ranges),
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_qanbari_kauai_ehh, 4))


## CLR


sim_bermuda_kauai_clr <- simulate_overlaps(width(sweep_ranges$CLR_bermuda) - 1,
                                           width(sweep_ranges$CLR_kauaii) - 1,
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_bermuda_kauai_clr, 44))


sim_qanbari_bermuda_clr <- simulate_overlaps(width(sweep_ranges$CLR_bermuda) - 1,
                                             width(qanbari_ranges),
                                             autosomal_size,
                                             n_rep = 1000)

print(analyse_sim(sim_qanbari_bermuda_clr, 10))


sim_qanbari_kauai_clr <- simulate_overlaps(width(sweep_ranges$CLR_kauai) - 1,
                                           width(qanbari_ranges),
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_qanbari_kauai_clr, 21))


## Tajima


sim_bermuda_kauai_taj <- simulate_overlaps(width(sweep_ranges$Tajima_Bermuda) - 1,
                                           width(sweep_ranges$Tajima_Hawaii) - 1,
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_bermuda_kauai_taj, 39))

sim_qanbari_bermuda_taj <- simulate_overlaps(width(sweep_ranges$Tajima_Bermuda) - 1,
                                             width(qanbari_ranges),
                                             autosomal_size,
                                             n_rep = 1000)

print(analyse_sim(sim_qanbari_bermuda_taj, 18))


sim_qanbari_kauai_taj <- simulate_overlaps(width(sweep_ranges$Tajima_Hawaii) - 1,
                                           width(qanbari_ranges),
                                           autosomal_size,
                                           n_rep = 1000)

print(analyse_sim(sim_qanbari_kauai_taj, 17))



## Consensus Bermuda / Hawaii


sim_consensus <- simulate_overlaps(width(consensus_ranges$bermuda) - 1,
                                   width(consensus_ranges$kauai) - 1,
                                   autosomal_size,
                                   n_rep = 1000)

print(analyse_sim(sim_consensus, 3))
