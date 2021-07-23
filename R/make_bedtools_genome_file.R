
library(Biostrings)


genome <- readDNAStringSet("galGal4_validated.fa")
chr_size <- data.frame(chr = names(genome), length = width(genome))
write.table(chr_size, "galGal4_chr_sizes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
