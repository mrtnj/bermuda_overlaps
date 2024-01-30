## Check overlaps of sweeps with structural variants detected with Delly

library(assertthat)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)


## Delly calls, raw

files <- system("ls data/merged/delly_merged_raw*.txt",
                intern = TRUE)
names(files) <- files

sv <- map_dfr(files,
              read_tsv,
              col_names = FALSE,
              na = ".")
colnames(sv)[1:5] <- c("chr", "start", "info", "qual", "filter")

sv <- extract(sv,
              info,
              into = "type",
              regex = "SVTYPE=(DEL|DUP|INS|INV)",
              remove = FALSE)

sv <- extract(sv,
              info,
              into = "end",
              regex = "END=([0-9]+)",
              remove = FALSE,
              convert = TRUE)


## Get allele frequencies

geno <- sv[, 8:ncol(sv)]

get_altcount <- function(geno) {

    assert_that(ncol(geno) == 21)
    
    geno_recoded <- map_dfc(geno, function(g) {
        case_when(g == "0/0" ~ 0,
                  g == "0/1" ~ 1,
                  g == "1/1" ~ 2,
                  g == "./." ~ NA_real_)
    })
    
    ac <- rowSums(geno_recoded, na.rm = TRUE)
    
    assert_that(length(ac) == nrow(geno))

    ac
}

ac <- get_altcount(geno)

## Allele frequency filtered

sv_filtered <- sv[ac > 1 & sv$filter == "PASS",]

print(table(sv_filtered$type))


sv_ranges <- makeGRangesFromDataFrame(sv_filtered,
                                      keep.extra.columns = TRUE)



## Bermuda sweeps

files <- system("ls data/Bermuda*_regions",
                intern = TRUE)

sweeps <- map(files, read_tsv)

sweep_ranges <- map(sweeps, makeGRangesFromDataFrame)

names(sweep_ranges) <- c("Hp", "Tajima's D", "Both Hp and Tajima's D")



overlaps_sv <- as.data.frame(subsetByOverlaps(sv_ranges, sweep_ranges[[3]]))
overlaps_sweep <- as.data.frame(subsetByOverlaps(sweep_ranges[[3]], sv_ranges))

overlaps_sv$altcount <- get_altcount(overlaps_sv[,grepl("^X", colnames(overlaps_sv))])



## Graph

overlaps_sv <- mutate(group_by(overlaps_sv, seqnames),
                      xpos = 1 + 0.4 * 1:length(start))
overlaps_sweep <- mutate(group_by(overlaps_sweep, seqnames),
                          xpos = 2 + 0.4 * 1:length(start))
plot_overlap_sv <- ggplot() +
    geom_linerange(aes(x = xpos,
                       ymin = start/1e6,
                       ymax = end/1e6,
                       colour = type),
                   size = 2,
                   data = overlaps_sv) +
    geom_text(aes(label = paste(signif(width/1e3, 3), "kbp"),
                  x = xpos - 0.1,
                  y = (start + end)/2/1e6,
                  colour = type),
              hjust = "inward",
              data = overlaps_sv) +
    geom_text(aes(label = paste(signif(width/1e3, 3), "kbp"),
                  x = xpos - 0.1,
                  y = (start + end)/2/1e6),
              data = overlaps_sweep) +
    geom_linerange(aes(x = xpos,
                       ymin = start/1e6,
                       max = end/1e6),
                   size = 2,
                   data = overlaps_sweep) +
    facet_wrap(~ seqnames, scale = "free_x") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    coord_flip() +
    ylab("Position (Mbp)") +
    xlab("") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ggtitle("Structural variant calls overlapping sweep regions")


pdf("figures/sweep_sv.pdf")
print(plot_overlap_sv)
dev.off()
