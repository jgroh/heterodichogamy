library(ggplot2)
library(data.table)
library(seqinr)
library(ape)

# ----- TPPD
gff <- read.gff("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/BNU/JRE.chrGene.gff3")
setDT(gff)
TPP <- gff[grepl('JRE010410.1', x = attributes, ignore.case = T) & type != 'mRNA', .(type, start, end)]

#hJregCha_TPPseq <- read.fasta("~/workspace/heterodichogamy/H_locus_structure/TPPD/CDS/h_Jreg_Cha_CDS.fasta")
#ref_cds_seq <- as.character(hJregCha_TPPseq[[1]])

# ----- regia SNP positions ----- 
leg <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2BNU_HJ1-0_TPPD_HJ1SNP_phased.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2BNU_HJ1-0_TPPD_HJ1SNP_phased.impute.hap")
setnames(haps0, paste0("V", 1:226), paste0(rep(seq(1:113), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2BNU_HJ1-0_TPPD_HJ1SNP_phased.impute.hap.indv", header = F, col.names = c("run"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]

# ----- genotypes, assigned by coverage -----
sra_geno <- fread("~/workspace/heterodichogamy/regia/sra_samples_geno.txt")
sra_geno[, ID := NULL]
sra_geno[, loc := NULL]

d2 <- merge(sra_geno, d1[run %in% sra_geno$run], by = 'run', all = T)

# ------ read phenotypes for founders -----
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("ID", "phenotype", "species"))

pheno[ID != 'JG0026' & phenotype == 'protogynous' & species == 'regia', genotype := 'Hh']
pheno[ID == 'JG0026', genotype := 'HH']
pheno[phenotype=='protandrous', genotype := 'hh']

reg_geno <- pheno[species == 'regia', .(run = ID, genotype)]

d3 <- merge(reg_geno, d1)
# d3[, length(unique(run))] # why is one missing?


# ------ final data set ------
haps <- rbind(d2, d3)
haps[, length(unique(run))]



# ----- Assign haplotype identities ----

haps[pos == min(pos) & allele == 0, haplotype := 'H']
haps[pos == min(pos) & allele == 1, haplotype := 'h']

# ----- plot haplotypes

haps[, SNPrank := seq(1, .N), by = .(run, hap.id)]

setkey(haps, haplotype)

haps[, hapRank := seq(1, .N), by = .(pos)]


# look at haplotypes across HJ1-0 and TPPD
ggplot(haps[pos < 30784255 & pos < 30784440 ], 
       aes(x = SNPrank, y = hapRank)) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(
    #axis.text = element_blank(),
        legend.position = 'none')


# across just TPPD, potentially no phasing errors?
ggplot(haps[pos > 30784032 & pos < 30786419], 
       aes(x = SNPrank, y = hapRank)) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')




# ========= Align with gene =====

# across just TPPD, potentially no phasing errors
hapsTPPD <- haps[pos > 30784032 & pos < 30786419]
hapsTPPD[, SNPrankTPPD := seq(1, .N), by = .(run, hap.id)]
hapsTPPD[, hapRankTPPD := seq(1, .N), by = .(pos)]

nSNP <- hapsTPPD[, max(SNPrankTPPD,na.rm=T)]

ggplot(hapsTPPD, 
       aes(x = SNPrankTPPD, y = hapRank)) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  scale_x_continuous(breaks = 1:115, labels = NULL, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(aspect.ratio = 0.3,
        axis.text = element_blank(),
        legend.position = 'none', 
        axis.line=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.05, 'cm')
  )

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
#lines(x = c(0, 1), y = c(0, 0))
lines(x = c(0, 1), y = c(.3, .3))

pos_scld <- unique(haps[pos > 30784032 & pos < 30786419, pos])
pos_scld <- (pos_scld - 30784032)/(max(TPP[, end]) - 30784032)

for(i in 1:nSNP){
  lines(x = c(i/nSNP, pos_scld[i]), y= c(0.5, 0.3), lwd = 0.5)
}

TPP[, start_scld := (start - min(start))/(max(end) - min(start))]
TPP[, end_scld := (end - min(start))/(max(end) - min(start))]

for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'exon']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.25, 0.3, 0.3, 0.25), col = 'gray')
  }
}
for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'CDS']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.25, 0.3, 0.3, 0.25), col = 'goldenrod3')
  }
}


# ===== Tally Polymorphisms and fixed differences =====

# restrict to CDS regions
cds_regions <- TPP[type == 'CDS']
setkey(cds_regions, start)

# get list of SNPs that fall within cds
cds_snp_pos <- NULL
for(p in unique(hapsTPPD$pos)){
  if(any(cds_regions[, p >= start & p <= end])){
    cds_snp_pos <- c(cds_snp_pos, p)
  }
}

cdsHaps <- hapsTPPD[pos %in% cds_snp_pos]

# get allele counts by haplotype at each site
cds_SNP_counts <- cdsHaps[, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
cds_SNP_counts <- dcast(cds_SNP_counts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')

# **reverse complement these to be on the correct strand**
cds_SNP_counts[, c('allele0', 'allele1') := lapply(.SD, comp), .SDcols = c('allele0', 'allele1')]

# how many of each haplotype?
hapsTPPD[pos == min(pos), table(haplotype)]

# get CDS sequence

CDS <- cds_regions[, seq(start,end), by = start_scld][, .(pos = rev(V1))]

CDS[, seq_rank := seq(1, .N)]
CDS[, refseq := ref_cds_seq]

cds_SNP_counts <- merge(CDS[, .(pos, seq_rank)], cds_SNP_counts)

setkey(cds_SNP_counts, seq_rank)
plot(cds_SNP_counts$seq_rank ~ cds_SNP_counts$pos)

cds_SNP_counts[seq_rank %% 3 == 1, codon_position := 1]
cds_SNP_counts[seq_rank %% 3 == 2, codon_position := 2]
cds_SNP_counts[seq_rank %% 3 == 0, codon_position := 3]


# ----- big function to classify as synonymous vs. nonsynonymous -----

# lookup is the following table
# pos rank refseq
# 1: 31887006    1      a
# 2: 31887005    2      t
# 3: 31887004    3      g

# the rank variable should access the same position in the CDS. 
# also, basepairs should be according to the transcript sequence (may be opposite strand from reference sequence)
# CDS should be keyed by rank

syn_vs_nonsyn <- function(rnk, nuc, codon_position){
  
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

cds_SNP_counts[, type := syn_vs_nonsyn(rnk = seq_rank, nuc = allele1, codon_position = codon_position), by = pos]


# make sure you double check final numbers when assigning these classes
cds_SNP_counts[H == 54 & h == 0, fixed_poly := 'fixed']
cds_SNP_counts[h > 0 & H < 54, fixed_poly := 'polymorphic_both']
cds_SNP_counts[h > 0 & H == 54, fixed_poly := 'polymorphic_h']
cds_SNP_counts[H < 54 & h == 0, fixed_poly := 'polymorphic_H']


fxd <-  ggplot(cds_SNP_counts[fixed_poly == 'fixed'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 7)) + 
  labs(title = 'Fixed', x = '', y = 'Count')  + 
  theme(aspect.ratio = 1)

polymorphic_H <- ggplot(cds_SNP_counts[fixed_poly == 'polymorphic_H'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 7)) + 
  labs(title = 'Polymorphic in H', x = '', y = '') + 
  theme(aspect.ratio = 1)

polymorphic_h <- ggplot(cds_SNP_counts[fixed_poly == 'polymorphic_h'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  labs(title = 'Polymorphic in h', x = '', y = '') + 
  scale_y_continuous(limits = c(0, 7)) + 
  theme(aspect.ratio = 1)

grid.arrange(fxd,polymorphic_H, polymorphic_h, ncol = 3)

divN <- cds_SNP_counts[fixed_poly == 'fixed' & type == 'N', .N]
divS <- cds_SNP_counts[fixed_poly == 'fixed' & type == 'S', .N]
polyN <- cds_SNP_counts[fixed_poly == 'polymorphic_H' & type == 'N', .N]
polyS <- cds_SNP_counts[fixed_poly == 'polymorphic_H' & type == 'S', .N]


z <- matrix(c(polyS, divS, polyN, divN), nrow = 2)

fisher.test(z)

# option to write file out, then samtools faidx the reference to get the cds sequence
# (had previously done, just read at top of script)
# fwrite(data.table(CDS[, paste0("NC_049911.1:", pos, "-", pos)]), 
#       file= "~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/TPPD_CDS_loc.txt",
#       quote = F, row.names = F, col.names = F)

