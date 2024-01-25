library(ape)


#x <- read.dna("~/workspace/heterodichogamy/01_GC_locus_structure/dnds_alignments/XBAT33-like/Lakota1_Jreg.fa", format = 'fasta')


x <- read.dna("~/workspace/heterodichogamy/01_GC_locus_structure/dnds_alignments/19_BAG-family-molecular-chaperone-regulator-1-like/Ccat_Jreg.fasta", format = 'fasta')


z <- dnds(x, code = 1, codonstart = 1, quiet = FALSE, detail = TRUE, return.categories = F)

#z
# number of degeneracy categories
n <- dnds(x, code = 1, codonstart = 1, quiet = FALSE, details = TRUE, return.categories = TRUE)

#length(which(n[1,] == 0))
#length(which(n[1,] == 2))
length(which(n[1,] == 4))



# ====== Date Estimate =====
JC <- function(p){
  -(3/4)*log(1 - (4/3)*p)
}

# G-loc divergence estimates
Carya_Gloc_div <- mean(c(0.059401884,0.046910755))
Gloc2regia_div <- mean(c(0.066722269,0.068631579,0.062068966,0.059390048))

S_G <- JC(Carya_Gloc_div)
S_reg <- JC(Gloc2regia_div)

(S_G/S_reg)*58
(S_G/S_reg)*72


