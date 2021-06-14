
## Overlap Bermuda sweeps with a lifted version of the Groenen et al. 2009
## linkage map, in order to understand whether there is an association between
## Bermuda sweeps and low recombination rate

library(dplyr)
library(GenomicRanges)
library(ggbeeswarm)
library(ggplot2)
library(readr)
library(purrr)


## Read sweeps

files <- system("ls data/Bermuda*_regions",
                intern = TRUE)

sweeps <- map(files, read_tsv)

sweep_ranges <- map(sweeps, makeGRangesFromDataFrame)

names(sweep_ranges) <- c("Hp", "Tajima's D", "Both Hp and Tajima's D")


## Read recombination regions

rec <- read_tsv("data/groenen2009_windows_500kbp_galgal4.txt")

rec_ranges <- makeGRangesFromDataFrame(rec, keep.extra.columns = TRUE)


get_sweep_rec <- function(ranges) {

    n_ranges <- length(ranges)
    
    sweep_rec <- numeric(n_ranges)
    
    for (sweep_ix in 1:n_ranges) {
        
        rec_overlaps <- subsetByOverlaps(rec_ranges,
                                         ranges[sweep_ix])
        
        sweep_rec[sweep_ix] <- mean(mcols(rec_overlaps)$average_rec)
        
    }
    
    data.frame(as.data.frame(ranges),
               rec = sweep_rec)
}


sweep_rec <- map_dfr(sweep_ranges,
                     get_sweep_rec,
                     .id = "set")

background <- as.data.frame(rec_ranges)[, 1:6]
colnames(background)[6] <- "rec"
background$set <- "Whole genome"


medians <- summarise(group_by(sweep_rec, set),
                     median = median(rec, na.rm = TRUE))


plot_rec_sweep <- ggplot() +
    geom_quasirandom(aes(x = set,
                         y = rec),
                     data = rbind(sweep_rec, background)) +
    geom_hline(yintercept = median(rec$average_rec, na.rm = TRUE),
               colour = "blue") +
    geom_point(aes(x = set,
                   y = median),
               colour = "red",
               data = medians) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ylab("Recombination rate (cM/Mb)") +
    xlab("") +
    ggtitle("Recombination rate around sweep regions")
    

wilcox.test(x = filter(sweep_rec, set == "Hp")$rec,
            y = rec$average_rec)

wilcox.test(x = filter(sweep_rec, set == "Tajima's D")$rec,
            y = rec$average_rec)

wilcox.test(x = filter(sweep_rec, set == "Both Hp and Tajima's D")$rec,
            y = rec$average_rec)
