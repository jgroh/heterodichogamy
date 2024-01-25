library(data.table)

#########################
# Currently this script uses the annotation from the BNU assembly, as reads were aligned to this to avoid copy number errors in phasing. 
# However the annotation for this gene appears to be only partial, so will keep original analysis. 



# ----- SNP positions ----- 
# BNU
#leg <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2BNU_chr7_30.75-30.8Mb_phased.impute.legend")

# Chandler
leg <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_TPPD1.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_TPPD1.impute.hap")
colnames(haps0)
n <- ncol(haps0)
setnames(haps0, paste0("V", 1:n), paste0(rep(seq(1:(n/2)), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_TPPD1.impute.hap.indv", header = F, col.names = c("ID"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]

d1


# ---- merge with genotype (assigned a priori or by coverage) -----
geno <- fread("~/workspace/heterodichogamy/regia/sra_samples_geno.txt")

pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("ID", "phenotype", 'species'), header = F)
pheno[species == 'regia' & phenotype == 'protandrous', genotype := 'gg']
pheno[species == 'regia' & phenotype == 'protogynous' & ID != 'JG0026', genotype := 'Gg']
pheno[species == 'regia' & phenotype == 'protogynous' & ID == 'JG0026', genotype := 'GG']

regia_genotypes <- rbind(geno[, .(ID=run, genotype)], pheno[species == 'regia', .(ID, genotype)])

haps <- merge( regia_genotypes, d1, by = 'ID', all = T)


# ----- visualize distribution of nonreference alleles -----
# (expect to be bimodal)
# BNU
#a <- haps[pos > 30784032 & pos < 30786419 & hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(ID,genotype)]
#b <- haps[pos > 30784032 & pos < 30786419 & hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(ID,genotype)]

# Chandler
a <- haps[pos > 31884475 & pos < 31887006 & hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(ID,genotype)]
b <- haps[pos > 31884475 & pos < 31887006 & hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(ID,genotype)]

nonref_alleles <- merge(a,b)
nonref_allhaps <- melt(nonref_alleles, id.vars = c("ID", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)


# ----- Assign haplotype ID based on fake SNP that codes for SV -----
haps[, hapRank := seq(1, .N), by = .(pos)]

#G.IDs <- haps[pos == 30780000 & allele == 0, hapRank]
#g.IDs <- haps[pos == 30780000 & allele == 1, hapRank]

g.IDs <- haps[pos == 31884475 & allele == 0, hapRank]
G.IDs <- haps[pos == 31884475 & allele == 1, hapRank]
haps[hapRank %in% G.IDs, hap := 'G']
haps[hapRank %in% g.IDs, hap := 'g']



# ----- visualize haplotypes -----
haps[, SNPrank := seq(1, .N), by = .(ID, hap.id)]


# BNU NDR1 coordinates 30759890-30760855

ggplot(haps[ pos > 31884475 & pos < 31887006], 
       aes(x = SNPrank, y = hapRank)) + 
  facet_wrap(~hap) +
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none') 

# which individual has haplotype that looks incorrectly phased?
haps[hap == 'g' & pos > 31884475 & pos < 31887006, sum(allele), by = ID][V1 > 10]
# This is from a heterozygote, looks like phase switches in that individual. Just discard. 

# subset haps to remove this individual and also only to regions within the TPPD exons
haps <- haps[ID != 'SRR14430327' & pos > 31884475 & pos < 31887006]

# look just at heterozygotes 
ggplot(haps[genotype == 'Gg' & 
              #pos > 30759890 & pos < 30786419],  #BNU
              pos > 31884475 & pos < 31887006],
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~ID) +
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')  



haps[genotype == 'gg',]



# =========== MK TEST ==================
# ====== Sequence Analysis =====
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


# coding sequence coordinates
#coords01 <- fread('~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/BNU/Jregia_BNU_TPPD-1.gff', sep = '\t')
coords01 <- fread("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/Jregia_TPPD1_CDS.gff", sep = '\t')
coords01 <- coords01[, c(4,5,7)]
setnames(coords01, c("start", "end", 'strand'))
setkey(coords01, "start")

#seq01 <- read.fasta('~/workspace/heterodichogamy/01_G_locus_structure/TPPD/CDS/Jregia_BNU_TPPD-1.fasta')
seq01 <- read.fasta('~/workspace/heterodichogamy/01_G_locus_structure/TPPD/CDS/h_Jreg_Cha_CDS.fasta')
seq01 <- as.character(seq01[[1]])
#length(seq01) 

# combine positions and CDS sequence
if(coords01[1, strand == '-']){
  cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = rev(V1))], refseq=seq01)
} else{
  cds01 <- data.table(coords01[, seq(start,end), by = start][, .(pos = V1)], refseq=seq01)
}

cds01


# assign codon position
cds01[, seq_rank := seq(1, .N)]
cds01[seq_rank %% 3 == 1, codon_position := 1]
cds01[seq_rank %% 3 == 2, codon_position := 2]
cds01[seq_rank %% 3 == 0, codon_position := 3]

# get allele counts by haplotype at each site
cnts01 <- haps[pos >= coords01[, min(start)] & pos <= coords01[, max(end)], .(count = sum(allele)), by = .(pos, hap, allele0, allele1)]
cnts01 <- dcast(cnts01, pos + allele0 + allele1 ~ hap, value.var = 'count')

# *conditional reverse complement*
if(coords01[1, strand == '-']){
  cnts01[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]
}

final01 <- merge(cnts01, cds01, by = 'pos')
setkey(final01, seq_rank)


# classify
vals <- vector()
for(i in 1:nrow(final01)){
  vals[i] <- syn_vs_nonsyn(CDS_dt = cds01, rnk = final01$seq_rank[i], nuc = final01$allele1[i], codon_position = final01$codon_position[i])
}

final01[, type := vals]

# classify polymorphisms
# how many of each morph/ genotype/ hap?
haps[SNPrank == 2,
     table(hap)]

haps[SNPrank == 2 & hap.id == 1,
     table(genotype)]

# chi squared test for G-locus alleles
HWChisq(X = c(AA=61,AB=50,BB=2))




# ----- classify polymorphisms -----
final01[(G >= 52 & g <= 1) | (G <= 1 & g >= 171), fixed_poly := 'fixed']
final01[G > 1 & G < 52 & g > 1 & g < 171, fixed_poly := 'polymorphic_both']
final01[(G > 1 & G < 52) & (g <= 1 | g >= 171), fixed_poly := 'polymorphic_G']
final01[(G <= 1 | G >= 52) & (g > 1 & g < 171), fixed_poly := 'polymorphic_g']


final01[, fixed_poly_MK := fixed_poly]
final01[fixed_poly == 'polymorphic_both', fixed_poly_MK := 'polymorphic_g']

MK_tbl <- table(final01$fixed_poly_MK, final01$type)
MK_tbl


# ====================
maf <- '~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments/CillPaw/chr11.maf'

# read alignment. Use sep value that's not in file to force fread to behave like readLines
mafLines <- fread(maf, sep = "?", skip = 1,
                  header = F)

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
mafDT <- mafDT[, .(alignmentID, chr, start0 = as.numeric(start0), length = as.numeric(length), seq)]
mafDT[, start1 := start0 + 1]
# mafDT[, c(1:4)]


# function that takes alignment block and returns data table with aligned seqs in long format. First row in alignment block is ref, 2nd is qry.
tr <- function(DT, by){
  outDT <- data.table(alignmentID = by,
                      refChr = DT[1, chr], 
                      qryChr = DT[2, chr], 
                      refStart = DT[1, start0], 
                      qryStart = DT[2, start0],
                      refLength = DT[1, length], 
                      qryLength = DT[2, length],
                      refSeq = unlist(strsplit(DT[1, seq], "")),
                      qrySeq = unlist(strsplit(DT[2, seq], ""))
  )
  
  # ignore sites in the alignment that are absent in the reference so that coordinates can be earily set for reference in next step
  outDT <- outDT[refSeq != "-"]
  
  # set coordinates of reference. + 1 is because MAF start position is zero-based
  outDT[, refPos := seq_len(.N)]
  outDT[, refPos := refPos + refStart]
  outDT
}

seqDT <- mafDT[, tr(DT = .SD, by = .BY), by = alignmentID]


TPPD1 <- seqDT[refPos >= 31884478 & refPos <= 31887006]
TPPD1

# -------- combine to polarize SNPs -----
polarized <- merge(final01, TPPD1[, .(pos = refPos, refSeq, qrySeq)], by = 'pos')
polarized

# only do if on reverse strand
polarized[, c('refSeq', 'qrySeq') := lapply(.SD, comp), .SDcols = c('refSeq', 'qrySeq')]
polarized
polarized <- polarized[fixed_poly == 'fixed']

# These are the SNPs where there's a fixed difference in TPP coding sequence between J. regia G-locus haplotypes
# and the reference (J. regia Chandler) matches the outgroup (pecan). So these are derived SNPs in the Juglans G lineage
polarized[qrySeq == allele0]

# These are SNPs where the outgroup matches the 'alternate' allele (G haplotypes), so the SNPs fixed in the Juglans g lineage
polarized[qrySeq == allele1]


# is there a multiple hit
p1 <- polarized[qrySeq == allele0, pos]
p2 <- polarized[qrySeq == allele1, pos]
polarized[!pos %in% c(p1, p2)]

