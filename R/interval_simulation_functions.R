library(IRanges)

## Function to place a random interval

random_interval <- function(interval_length,
                            genome_size) {
    
    start <- sample(1:(genome_size - interval_length), 1)
    
    IRanges(start = start,
            end = start + interval_length - 1)
}

## Function to place a set of random intervals 

random_interval_set <- function(interval_lengths,
                                genome_length) {
    
    sim_intervals <- lapply(interval_lengths,
                            random_interval,
                            genome_size = autosomal_size)
    
    Reduce(c, sim_intervals)
    
}



simulate_overlaps <- function(widths_set1,
                              widths_set2,
                              autosomal_length,
                              n_rep) {
    
    sim_set1 <- replicate(n_rep,
                          random_interval_set(widths_set1,
                                              autosomal_length),
                          simplify = FALSE)
    
    sim_set2 <- replicate(n_rep,
                          random_interval_set(widths_set2,
                                              autosomal_length),
                          simplify = FALSE)
    
    
    ## Count overlaps between simulated Bermuda and Kauai, Bermuda and domestic
    
    count_set_overlaps <- function(x, y) length(subsetByOverlaps(x, y))
    
    sim_overlaps <- mapply(count_set_overlaps, sim_set1, sim_set2)
    
    sim_overlaps
}