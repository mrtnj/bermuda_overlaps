## Look at coverage data from sweep regions

library(dplyr)
library(GenomicRanges)
library(ggbeeswarm)
library(ggplot2)
library(readr)
library(purrr)

## Coverage data

files <- dir("outputs/coverage")

sample_ids <- sub(files,
                  pattern = "_coverage.txt",
                  replacement = "")

files <- paste("outputs/coverage/", files, sep = "")

names(files) <- sample_ids

coverage <- map_dfr(files,
                    col_names = FALSE,
                    read_tsv,
                    .id = "sample_id")

colnames(coverage) <- c("sample_id", "chr", "start", "end",
                        "overlapping_reads", "bases_covered",
                        "window_length", "fraction_covered")

coverage_ranges <- GRanges(seqnames = coverage$chr,
                           ranges = IRanges(coverage$start, coverage$end),
                           mcols = coverage)


## Bermuda sweeps

files <- system("ls data/Bermuda*_regions",
                intern = TRUE)

sweeps <- map(files, read_tsv)

sweep_ranges <- map(sweeps, makeGRangesFromDataFrame)

names(sweep_ranges) <- c("Hp", "Tajima's D", "Both Hp and Tajima's D")



## Average coverage of sweeps

average_coverage <- summarise(group_by(coverage, chr, start, end),
                              average_coverage = mean(overlapping_reads))

average_coverage_ranges <- makeGRangesFromDataFrame(average_coverage, 
                                                    keep.extra.columns = TRUE)


average_coverage_sweeps <- map_dfr(sweep_ranges,
                                   function(sr) {
                                       overlap <- subsetByOverlaps(average_coverage_ranges,
                                                                   sr)
                                       
                                       overlap_df <- as.data.frame(overlap)
                                       
                                       overlap_df <- overlap_df[, c("seqnames", "start", "end", "average_coverage")]
                                       colnames(overlap_df)[1] <- "chr"
                                       
                                       overlap_df
                                   },
                                   .id = "set")

background <- transform(average_coverage,
                        set = "Whole genome")
background <- background[, c("set", "chr", "start", "end", "average_coverage")]


medians <- summarise(group_by(average_coverage_sweeps, set),
                     median = median(log10(average_coverage), na.rm = TRUE))

plot_coverage_sweep <- ggplot() +
    geom_quasirandom(aes(x = set,
                         y = log10(average_coverage)),
                     data = rbind(average_coverage_sweeps, background)) +
    geom_hline(yintercept = median(log10(background$average_coverage), na.rm = TRUE),
               colour = "blue") +
    geom_point(aes(x = set,
                   y = median),
               colour = "red",
               data = medians) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ylab("Logarithm of number of reads (10 kbp bins)") +
    xlab("") +
    ggtitle("Coverage around sweep regions")

pdf("figures/sweep_coverage.pdf",
    height = 3.5)
print(plot_coverage_sweep)
dev.off()


wilcox.test(x = filter(average_coverage_sweeps, set == "Both Hp and Tajima's D")$average_coverage,
            y = background$average_coverage)
