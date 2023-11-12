library(ape)

x <- read.dna("~/workspace/heterodichogamy/H_locus_structure/TPPD/transcripts/aligned_transcripts.fasta", format = 'fasta')


# value of dNdS
dnds(x, code = 1, codonstart = 1, quiet = FALSE, details = TRUE, return.categories = FALSE)

# number of degeneracy categories
n <- dnds(x, code = 1, codonstart = 1, quiet = FALSE, details = TRUE, return.categories = TRUE)

length(which(n[1,] == 0))
length(which(n[1,] == 2))
length(which(n[1,] == 4))


#4/179 

4/179 
