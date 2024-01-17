library(ggplot2)
library(data.table)
library(seqinr)
library(ape)


# ========== Format alignment result between pecan 'Pawnee' and J. regia V2 assemblies ======

maf <- '~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments/CillLak1/Carya_chr4.maf'

# read alignment. Use sep value that's not in file to force fread to behave like readLines
mafLines <- fread(maf, sep = "?", skip = 1,
                  header = F)


# ===== format maf alignment into data.table =====

# subset to alignment blocks
mafLines <- mafLines[grep("^s", V1)]

# set column names
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5
mafDT <- mafLines[, tstrsplit(V1, "\\s+")]
setnames(mafDT, c("s", "chr", "start0", "length", "strand", "srcLength", "seq"))


# add grouping variable for alignment id
mafDT[, alignmentID := rep(1:(.N / 2), each = 2)]
# mafDT[, c(1:6,8)]

# subset to focal columns
mafDT <- mafDT[, .(alignmentID, chr, start0 = as.numeric(start0), length = as.numeric(length), strand, seq)]
mafDT[, start1 := start0 + 1]
# mafDT[, c(1:4)]

# *swqp reference and query* so that pecan is the reference here
x <- c(2,1)
for(i in 1:((nrow(mafDT) - 1)/2)){ x <- c(x, c(2,1) + 2*i) }
mafDT <- mafDT[x]


mafDT[, c(1:5)]



# function that takes alignment block and returns data table with aligned seqs in long format. First row in alignment block is ref, 2nd is qry.
# definitely check this function for each specific use case.
# For example here I've swapped reference and query in maf after the alignment, code specific to this application .
tr <- function(DT, by){
  
  # deal with strand
  if(DT[1, strand == '-']){
    
    ref.Seq = rev(comp(unlist(strsplit(DT[1, seq], ""))))
    qry.Seq = rev(comp(unlist(strsplit(DT[2, seq], ""))))

    } else if(DT[1, strand == '+']){
    
      ref.Seq = unlist(strsplit(DT[1, seq], ""))
      qry.Seq = unlist(strsplit(DT[2, seq], ""))
      
    }
  
  
  ref.Start = DT[1, start1] 
  qry.Start = DT[2, start1] 
  
    
  outDT <- data.table(alignmentID = by,
                      refChr = DT[1, chr], 
                      qryChr = DT[2, chr], 
                      refStart = ref.Start, 
                      qryStart = qry.Start,
                      refLength = DT[1, length], 
                      qryLength = DT[2, length],
                      refSeq = ref.Seq,
                      qrySeq = qry.Seq
  )
  
  # ignore sites in the alignment that are absent in the reference so that coordinates can be earily set for reference in next step
  outDT <- outDT[refSeq != "-"]
  
  # set coordinates of reference. 
  
  
  outDT[, refPos := seq_len(.N)]
  outDT[, refPos := refPos + refStart - 1]
  outDT
  
}

seqDT <- mafDT[, tr(DT = .SD, by = .BY), by = alignmentID]


# ====== Tally Substitutions ======
# assign logical to substitutions
seqDT[, sub := ifelse(refSeq != qrySeq, T, F)]

# gaps in the query sequence are treated as NA 
# (we account for the total number of non-gapped bps in the calc of Dxy)
seqDT[qrySeq == "-", sub := NA]

# check total divergence 
# seqDT[, sum(sub, na.rm=T)/length(which(!is.na(sub)))]



# ======= Read phased haplotypes ======



# ----- SNP positions ----- 
leg <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_Lakota_v1_Gloc_CDS_phased.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_Lakota_v1_Gloc_CDS_phased.impute.hap")
colnames(haps0)
setnames(haps0, paste0("V", 1:60), paste0(rep(seq(1:30), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_Lakota_v1_Gloc_CDS_phased.impute.hap.indv", header = F, col.names = c("run"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]


# ---- merge with genotype (assigned by phenotype or coverage) -----
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")
haps <- merge(geno[, .(run =ID, genotype)], d1, by = 'run', all = T)
haps[genotype == 'hh', genotype := 'gg']
haps[genotype == 'Hh', genotype := 'Gg']
haps[genotype == 'HH', genotype := 'GG']


# ----- read Jregia ----- *replaced by long read alignment
# Jreg012 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Jreg_to_Lak1_invariant_sites.012")
# Jreg.pos <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Jreg_to_Lak1_invariant_sites.012.pos")
# Jreg.indv <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Jreg_to_Lak1_invariant_sites.012.indv", header = F)
# 
# setnames(Jreg012, c('ID', paste0("pos", Jreg.pos[, V2])))
# Jreg012[, ID := Jreg.indv[, 1]]
# 
# length(unique(Jreg.pos[, V2]))
# 
# 
# Jreg_geno <- melt(Jreg012, id.vars = 'ID')
# 
# Jreg_geno[, pos := rep(Jreg.pos[,V2], each = 3)]
# Jreg_geno[value == -1, value := NA]
# Jreg_allele <- Jreg_geno[, .(geno_sum = sum(value)), by = pos]
# 
# Jreg_allele[geno_sum == 0, Jreg_allele :=  'ref']
# Jreg_allele[geno_sum == 6, Jreg_allele :=  'alt']


# ----- visualize distribution of nonreference alleles -----
# (expect to be bimodal)
a <- haps[hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(run,genotype)]
b <- haps[hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(run,genotype)]
nonref_alleles <- merge(a,b)
nonref_allhaps <- melt(nonref_alleles, id.vars = c("run", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)

# just heterozygotes 
nonref_alleles[genotype == 'Gg', hap_more := max(nonref_alleles_hap1, nonref_alleles_hap2), by = .(run,genotype)]
nonref_alleles[genotype == 'Gg', hap_less := min(nonref_alleles_hap1, nonref_alleles_hap2), by = .(run,genotype)]

nonref <- melt(nonref_alleles[, .(run, hap_more, hap_less, genotype)], id.vars = c('run', 'genotype'), value.name = 'nonref', variable.name = 'haplotype')
ggplot(nonref[genotype == 'Gg'], aes(x = nonref)) + geom_histogram()
ggplot(nonref[genotype == 'Gg'], aes(x = haplotype, y = nonref, group = run)) + geom_line()


# ----- Assign haplotype ID based on fake SNP that codes for SV -----
haps[, hapRank := seq(1, .N), by = .(pos)]

G.IDs <- haps[pos == 6700000 & allele == 0, hapRank]
g.IDs <- haps[pos == 6700000 & allele == 1, hapRank]
haps[hapRank %in% G.IDs, hap := 'G']
haps[hapRank %in% g.IDs, hap := 'g']


# ----- visualize haplotypes -----
haps[, SNPrank := seq(1, .N), by = .(run, hap.id)]
haps[, hapRank := seq(1, .N), by = .(pos)]

ggplot(haps, 
       aes(x = SNPrank, y = hapRank)) + 
  facet_wrap(~hap) +
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')


# ====== Dxy =====
# 35463 total CDS base pairs
# index of all pairwise comparisons

# g haplotypes
# 
# gmat <- combn(g.IDs, 2)
# 
# g.pi.vec <- vector()
# for(i in 1:ncol(gmat)){
#   g.hap.pair <- dcast(data = haps[hapRank %in% c(gmat[,i]), .(pos, allele, hapRank)], 
#                       formula = pos ~ hapRank, value.var = 'allele')
#   setnames(g.hap.pair, c("pos", "hap1", "hap2"))
#   # calculate per bp pairwise differences
#   g.pi.vec[i] <- g.hap.pair[, sum(hap1!=hap2)/35463]
# }
# 
# g.pi <- mean(g.pi.vec) # 0.04161663
# g.pi
# 
# # ----- G haplotypes ----- 
# Gmat <- combn(G.IDs, 2)
# 
# G.pi.vec <- vector()
# for(i in 1:ncol(Gmat)){
#   G.hap.pair <- dcast(data = haps[hapRank %in% c(Gmat[,i]), .(pos, allele, hapRank)], 
#                       formula = pos ~ hapRank, value.var = 'allele')
#   setnames(G.hap.pair, c("pos", "hap1", "hap2"))
#   # calculate per bp pairwise differences
#   G.pi.vec[i] <- G.hap.pair[, sum(hap1!=hap2)/35463]
# }
# 
# G.pi <- mean(G.pi.vec) # 0.0068706
# G.pi
# 
# hist(G.pi.vec)
# hist(g.pi.vec)


# ----- permutation test ------

# # do for each ?????
#
# 
# ids.test <- sample(1:60)
# 
# G.ids.test <- 
#   
# haps[, hapRank := seq(1, .N), by = .(pos)]
# 
# haps[]
# sample(1:60)

# ========== Fixed Diffs =========


cnts <- haps[, .(cnt = sum(allele)), by = .(pos, hap, allele0, allele1)]
setkey(cnts, "pos")
cnts <- dcast(cnts, pos + allele0 + allele1 ~ hap, value.var = 'cnt')

cnts[, diff := abs(G-g)]
cnts[diff >= 39]

polarizedDiffs <- merge(cnts[diff >= 39], seqDT[, .(refSeq, qrySeq, pos = refPos)], by = 'pos')

# exclude dummy SNP
polarizedDiffs <- polarizedDiffs[pos != 6700000]

# double check that sequence indexing worked correctly
polarizedDiffs[refSeq != allele0]

polarizedDiffs[qrySeq == allele0, lineage := 'g']
polarizedDiffs[qrySeq == allele1, lineage := 'G']
polarizedDiffs[qrySeq != allele0 & qrySeq != allele1, lineage := 'ambiguous']
polarizedDiffs[, table(lineage)]


# tally for total substitutions
# 235 SNPs fixed in G lineage, 215 fixed in g lineage, 15 ambiguous ancestral state 


######################
# ====== MK TEST =====
######################


##############################################################
# ----- function to classify synonymous vs nonsynonymous -----
##############################################################


syn_vs_nonsyn <- function(CDS_dt, rnk, nuc, codon_position){
  CDS <- copy(CDS_dt)
  if(codon_position == 1){
    refaa <- translate(CDS[seq_rank %in% rnk:(rnk+2), refseq])
    qryaa <- translate(c(nuc, CDS[seq_rank %in% (rnk+1):(rnk+2), refseq]))
  }
  if(codon_position == 2){
    refaa <- translate(CDS[seq_rank %in% (rnk-1):(rnk+1), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-1), refseq], nuc, CDS[seq_rank %in% (rnk+1), refseq]) )
  }
  if(codon_position == 3){
    refaa <- translate(CDS[seq_rank %in% (rnk-2):(rnk), refseq])
    qryaa <- translate( c( CDS[seq_rank %in% (rnk-2):(rnk-1), refseq], nuc ) )
  }
  if(refaa == qryaa){
    return('S')
  } else if(refaa != qryaa){
    return('N')
  }
}

##################################
# ----- MK test function -----
# returns a table of nonsyn vs.
# syn fixed diffs and polymorphims
##################################

# forward strand test case
# fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/01_Protein_phosphatase_inhibitor_2-like/01.fasta'
# gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/01_Protein_phosphatase_inhibitor_2-like/01.gff'

# reverse strand test case
# fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/05_XBAT33-like/05.fasta'
# gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/05_XBAT33-like/05.gff'

MK_tbl <- function(fasta, gff){
  seq01 <- read.fasta(fasta)
  seq01 <- as.character(seq01[[1]])
  #length(seq01) 
  
  # coding sequence coordinates
  coords01 <- fread(gff, sep = '\t')
  coords01 <- coords01[, c(4,5,7)]
  setnames(coords01, c("start", "end", 'strand'))
  setkey(coords01, "start")
  
  
  # combine positions and CDS sequence
  if(coords01[1, strand == '-']){
    cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=seq01)
  } else{
    cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = V1)], refseq=seq01)
  }
  
  # assign codon position
  cds01[, seq_rank := seq(1, .N)]
  cds01[seq_rank %% 3 == 1, codon_position := 1]
  cds01[seq_rank %% 3 == 2, codon_position := 2]
  cds01[seq_rank %% 3 == 0, codon_position := 3]

  # get allele counts by haplotype at each site
  cnts01 <- haps[pos >= coords01[, min(start)] & pos <= coords01[, max(end)], .(count = sum(allele)), by = .(pos, hap, allele0, allele1)]
  cnts01 <- dcast(cnts01, pos + allele0 + allele1 ~ hap, value.var = 'count')
  
  final01 <- merge(cnts01, cds01, by = 'pos')
  setkey(final01, seq_rank)
  
  
  final01 <- merge(final01, seqDT[, .(refSeq, qrySeq, pos = refPos)], by = 'pos')
  
  # *conditional reverse complement*
  if(coords01[1, strand == '-']){
    final01[, c('allele0', 'allele1', 'refSeq', 'qrySeq') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1', 'refSeq', 'qrySeq')]
  }
  
  
  # classify
  vals <- vector()
  for(i in 1:nrow(final01)){
    vals[i] <- syn_vs_nonsyn(CDS_dt = cds01, rnk = final01$seq_rank[i], nuc = final01$allele1[i], codon_position = final01$codon_position[i])
  }
  
  final01[, type := vals]
  
  final01[abs(G-g) >= 39, fixed_poly := 'fixed']
  final01[G > 1 & G < 18 & g > 1 & g < 40, fixed_poly := 'polymorphic_both']
  final01[(G > 1 & G < 18) & (g <= 1 | g >= 40), fixed_poly := 'polymorphic_G']
  final01[(G <= 1 | G >= 18) & (g > 1 & g < 40), fixed_poly := 'polymorphic_g']
  
  
  final01[, fixed_poly_MK := fixed_poly]
  final01[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_g']
  
  
  MK_tbl <- table(final01$fixed_poly_MK, final01$type)
  #return(MK_tbl01)
  
  # which lineage is the variant derived in
  final01[qrySeq == allele0, lineage := 'g']
  final01[qrySeq == allele1, lineage := 'G']
  final01[qrySeq != allele0 & qrySeq != allele1, lineage := 'ambiguous']
  
  
  polarized_tbl <- final01[fixed_poly == 'fixed' & type == 'N', table(lineage)]
  
  return(list(fixed_poly_nonsyn_syn = MK_tbl, fixed_in_x_lineage = polarized_tbl))
}


cnt <- 0
tot <- 20

# ----- 01 -----
z01 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/01_Protein_phosphatase_inhibitor_2-like/01.fasta', 
        gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/01_Protein_phosphatase_inhibitor_2-like/01.gff')
z01

tot <- tot - 1
fisher.test(z01[[1]])




# ----- 02 -----
z02 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/02_LOC122306891_uncharacterized/02.fasta',
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/02_LOC122306891_uncharacterized/02.gff')
z02
cnt <- cnt + 1
fisher.test(z02[[1]])



# ----- 03 -----
z03 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/03_EMS1-like/03.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/03_EMS1-like/03.gff')
z03
fisher.test(z03[[1]])




# ----- 04 -----
z04 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/04_CEN-like-protein1/04.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/04_CEN-like-protein1/04.gff')
z04
fisher.test(z04[[1]])



# ----- 05 -----
z05 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/05_XBAT33-like/05.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/05_XBAT33-like/05.gff')
z05
cnt <- cnt + 1
fisher.test(z05[[1]])



# ----- 06 -----
z06 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/06_ROOT-PRIMORDIUM-DEFECTIVE-1/06.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/06_ROOT-PRIMORDIUM-DEFECTIVE-1/06.gff')
z06
cnt <- cnt + 1
fisher.test(z06[[1]])



# ----- 07 -----
z07 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/07_RHOMBOID-like-protein2/07.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/07_RHOMBOID-like-protein2/07.gff')
z07
fisher.test(z07[[1]])



# ----- 08 -----
z08 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/08_LOC122306943/08.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/08_LOC122306943/08.gff')
z08
fisher.test(z08[[1]])



# ----- 09 -----
z09 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/09_non-specific-lipid-transfer-protein14/09.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/09_non-specific-lipid-transfer-protein14/09.gff')
z09
cnt <- cnt + 1
fisher.test(z09[[1]])



# ----- 10 -----
z10 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/10_IQ-DOMAIN-23-like/10.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/10_IQ-DOMAIN-23-like/10.gff')
z10
cnt <- cnt + 1
fisher.test(z10[[1]])



# ----- 11 -----
z11 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/11_stamen-specific-protein-FIL1-like/11.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/11_stamen-specific-protein-FIL1-like/11.gff')
z11
cnt <- cnt + 1
fisher.test(z11[[1]])



# ----- 12 -----
z12 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/12_PMT16/12.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/12_PMT16/12.gff')
z12
cnt <- cnt + 1
fisher.test(z12[[1]])



# ----- 13 -----
z13 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/13_DNA-directed-RNA-polymerase-II-IV-V-subunit-6A-like/13.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/13_DNA-directed-RNA-polymerase-II-IV-V-subunit-6A-like/13.gff')
z13
fisher.test(z13[[1]])



# ----- 14 -----
z14 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/14_F-box-protein-PP2-A12-like/14.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/14_F-box-protein-PP2-A12-like/14.gff')
z14
cnt <- cnt + 1
fisher.test(z14[[1]])



# ----- 15 -----
z15 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/15_F-box-protein-At4g00755-like/15.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/15_F-box-protein-At4g00755-like/15.gff')
z15
#cnt <- cnt + 1
fisher.test(z15[[1]])



# ----- 16 -----
z16 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/16_SLK2/16.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/16_SLK2/16.gff')
z16
fisher.test(z16[[1]])



# ----- 17 -----
z17 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/17_protein-COFACTOR-ASSEMBLY-OF-COMPLEX-C-SUBUNIT-B-CCB2-chloroplastic-like/17.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/17_protein-COFACTOR-ASSEMBLY-OF-COMPLEX-C-SUBUNIT-B-CCB2-chloroplastic-like/17.gff')
z17
cnt <- cnt + 1
fisher.test(z17[[1]])



# ----- 18 -----
z18 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/18_dihydrodipicolinate-reductase-like-protein-CRR1-chloroplastic/18.fasta', 
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/18_dihydrodipicolinate-reductase-like-protein-CRR1-chloroplastic/18.gff')
z18
fisher.test(z18[[1]])



# ----- 19 -----
z19 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/19_BAG-family-molecular-chaperone-regulator-1-like/19.fasta',
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/19_BAG-family-molecular-chaperone-regulator-1-like/19.gff')
z19
cnt <- cnt + 1
fisher.test(z19[[1]])



# ----- 20 -----
z20 <- MK_tbl(fasta = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/20_protein-DETOXIFICATION-49-like/20.fasta',
              gff = '~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_genes/20_protein-DETOXIFICATION-49-like/20.gff')
z20
fisher.test(z20[[1]])


cnt



# how many nonsynonymous have fixed in g lineage
sum(
  z01[[2]]['g'], z02[[2]]['g'], z03[[2]]['g'], z04[[2]]['g'], z05[[2]]['g'],
  z06[[2]]['g'], z07[[2]]['g'], z08[[2]]['g'], z09[[2]]['g'], z10[[2]]['g'],
  z11[[2]]['g'], z12[[2]]['g'], z13[[2]]['g'], z14[[2]]['g'], z15[[2]]['g'], 
  z16[[2]]['g'], z17[[2]]['g'], z18[[2]]['g'], z19[[2]]['g'], z20[[2]]['g'],
  na.rm=T
)

# how many nonsynonymous have fixed in G lineage
sum(
  z01[[2]]['G'], z02[[2]]['G'], z03[[2]]['G'], z04[[2]]['G'], z05[[2]]['G'],
  z06[[2]]['G'], z07[[2]]['G'], z08[[2]]['G'], z09[[2]]['G'], z10[[2]]['G'],
  z11[[2]]['G'], z12[[2]]['G'], z13[[2]]['G'], z14[[2]]['G'], z15[[2]]['G'], 
  z16[[2]]['G'], z17[[2]]['G'], z18[[2]]['G'], z19[[2]]['G'], z20[[2]]['G'],
  na.rm=T
)


# tally for total substitutions
# 235 SNPs fixed in G lineage, 215 fixed in g lineage, 15 ambiguous ancestral state 

# tally for nonsynonymous substitutions
# 118 SNPs fixed in G lineage, 103 in g lineage


nG <- 118
sG <- 235-118

ng <- 103
sg <- 215-103

fisher.test(matrix(c(nG, ng, sG, sg), nrow=2))



