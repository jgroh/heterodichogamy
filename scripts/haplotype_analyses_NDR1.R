library(ggplot2)
library(data.table)
library(seqinr)
library(ape)

# ----- NDR1
gff <- read.gff("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/GCF_001411555.2_Walnut_2.0_genomic.gff")
setDT(gff)
NDR1 <- gff[grepl('XM_018957011.2', x = attributes, ignore.case = T) & type != 'mRNA', .(type, start, end)]

hJregCha_NDR1seq <- read.fasta("~/workspace/heterodichogamy/H_locus_structure/NDR1/h_Jreg_NDR1_CDS.fasta")
ref_cds_seq <- as.character(hJregCha_NDR1seq[[1]])

# ----- regia SNP positions ----- 
leg <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_NDR1_CDS_fakeSNPhJ2_phased.impute.legend")

# ID column in this table is both redundant and unecessary
leg[, ID := NULL]

# ----- haps ----- Each column is a biallelic SNP position. Columns are haplotypes.
haps0 <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_NDR1_CDS_fakeSNPhJ2_phased.impute.hap")
setnames(haps0, paste0("V", 1:226), paste0(rep(seq(1:113), each = 2), c("_1", "_2")))

d0 <- melt(cbind(leg, haps0), id.vars = c("pos", "allele0", "allele1"), value.name = 'allele')
d0[, c("indiv.id", "hap.id") := tstrsplit(variable, split = "_", fixed = T)]
d0[, variable := NULL]

# ----- individual IDs -----
indiv <- fread("~/workspace/heterodichogamy/regia/calls/Jregia_ALL2Cha_NDR1_CDS_fakeSNPhJ2_phased.impute.hap.indv", header = F, col.names = c("run"))
indiv[, indiv.id := as.character(seq(1, .N))]

d1 <- merge(d0, indiv, by = 'indiv.id', all = T)
d1[, indiv.id := NULL]

# ----- genotypes, assigned by coverage -----
sra_geno <- fread("~/workspace/heterodichogamy/regia/sra_samples_geno.txt")
sra_geno[, ID := NULL]
sra_geno[, loc := NULL]

d2 <- merge(sra_geno, d1[run %in% sra_geno$run], by = 'run', all = T)

# ------ read phenotypes for founders -----
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", header = F, col.names = c("ID", "phenotype", "species"))

pheno[ID != 'JG0026' & phenotype == 'protogynous' & species == 'regia', genotype := 'Hh']
pheno[ID == 'JG0026', genotype := 'HH']
pheno[phenotype=='protandrous', genotype := 'hh']

reg_geno <- pheno[species == 'regia', .(run = ID, genotype)]

d3 <- merge(reg_geno, d1)
# d3[, length(unique(run))] # why is one missing?


# ------ final data set ------
haps <- rbind(d2, d3)
haps[, length(unique(run))]

# ------ investigate error rate ------
# how many non-reference alleles does Chandler have? Should be zero. zero and zero, out of 64 SNPs.
haps[run=='Chandler', sum(allele), by = hap.id]
# leg[, length(unique(pos))]


# ----- investigate distribution of nonreference alleles among haplotypes

# we potentially expect to see difference in the number of nonref alleles between haplotypes
a <- haps[hap.id == 1, .(nonref_alleles_hap1 = sum(allele)), by = .(run,genotype)]
b <- haps[hap.id == 2, .(nonref_alleles_hap2 = sum(allele)), by = .(run,genotype)]
nonref_alleles <- merge(a,b)
nonref_allhaps <- melt(nonref_alleles, id.vars = c("run", 'genotype'), value.name = "Non ref alleles")

ggplot(nonref_allhaps, aes(x = `Non ref alleles`)) + geom_histogram() + 
  theme_classic() + 
  theme(aspect.ratio = 1)


# specifically look at heterozygotes
nonref_alleles[, hap_more := max(nonref_alleles_hap1, nonref_alleles_hap2), by = .(run,genotype)]
nonref_alleles[, hap_less := min(nonref_alleles_hap1, nonref_alleles_hap2), by = .(run,genotype)]

nonref <- melt(nonref_alleles[, .(run, hap_more, hap_less, genotype)], id.vars = c('run', 'genotype'), value.name = 'nonref', variable.name = 'haplotype')
ggplot(nonref[genotype == 'Hh'], aes(x = nonref)) + geom_histogram()
ggplot(nonref[genotype == 'Hh'], aes(x = haplotype, y = nonref, group = run)) + geom_line()


# ----- Assign haplotype identities ----
# assign based on whether a haplotype contains the hJ2 insertion (fake SNP)

tmp <- haps[pos == max(pos) & allele == 1, .(run, hap.id)]
tmp[, haplotype := 'H']

haps <- merge(haps, tmp, all.x = T)

# remaining ones are h haplotypes
haps[is.na(haplotype), haplotype := 'h']

# ----- plot haplotypes
haps[, SNPrank := seq(1, .N), by = .(run, hap.id)]

setkey(haps, haplotype, pos)

haps[, hapRank := seq(1, .N), by = .(pos)]

# look at haplotypes 
ggplot(haps[], 
       aes(x = SNPrank, y = hapRank)) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none')

# look at heterozygotes separately - 
ggplot(haps[genotype == 'Hh'], 
       aes(x = SNPrank, y = hap.id)) + 
  facet_wrap(~run) + 
  geom_tile(aes(fill = as.factor(allele)), width = 1, height = 1) + 
  scale_fill_manual(values = c('cornsilk', 'chocolate4')) + 
  theme_classic() + 
  labs(x = '', y = '', fill = 'Allele') + 
  theme(axis.text = element_blank(),
        legend.position = 'none', aspect.ratio = .1)




# ========= Align with gene =====

nSNP <- haps[, max(SNPrank,na.rm=T)]

ggplot(haps, 
       aes(x = SNPrank, y = hapRank)) + 
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



# ===== Tally Polymorphisms and fixed differences =====

# restrict to CDS regions
cds_regions <- NDR1[type == 'CDS']
setkey(cds_regions, start)

# get list of SNPs that fall within cds
cds_snp_pos <- NULL
for(p in unique(haps$pos)){
    if(any(cds_regions[, p >= start & p <= end])){
      cds_snp_pos <- c(cds_snp_pos, p)
    }
}

cdsHaps <- haps[pos %in% cds_snp_pos]

# get allele counts by haplotype at each site
cds_SNP_counts <- cdsHaps[, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
cds_SNP_counts <- dcast(cds_SNP_counts, pos + allele0 + allele1 ~ haplotype, value.var = 'count')


# get CDS sequence

CDS <- data.table(cds_regions[, seq(start,end)])[, .(pos = V1)]

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
  scale_y_continuous(limits = c(0, 10), breaks = 1:10) + 
  labs(title = 'Fixed', x = '', y = 'Count')  + 
  theme(aspect.ratio = 1)
  
polymorphic_H <- ggplot(cds_SNP_counts[fixed_poly == 'polymorphic_H'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 10), breaks = 1:10) + 
  labs(title = 'Polymorphic in H', x = '', y = '') + 
  theme(aspect.ratio = 1)

polymorphic_h <- ggplot(cds_SNP_counts[fixed_poly == 'polymorphic_h'], aes(x = type)) + 
  geom_bar(stat = 'count')  + 
  theme_classic() + 
  labs(title = 'Polymorphic in h', x = '', y = '') + 
  scale_y_continuous(limits = c(0, 10), breaks = 1:10) + 
  theme(aspect.ratio = 1)

grid.arrange(fxd,polymorphic_H, polymorphic_h, ncol = 3)

divN <- cds_SNP_counts[fixed_poly == 'fixed' & type == 'N', .N]
divS <- cds_SNP_counts[fixed_poly == 'fixed' & type == 'S', .N]
polyN <- cds_SNP_counts[fixed_poly == 'polymorphic_H' & type == 'N', .N]
polyS <- cds_SNP_counts[fixed_poly == 'polymorphic_H' & type == 'S', .N]


z <- matrix(c(polyS, divS, polyN, divN), nrow = 2)

fisher.test(z)


# ========= test whether nonsynonymous polymorphism is greater in H haplotypes
# ******* I don't think this approach is actually correct because it ignores the fact that these are two separate 'populations'
# strategy: permute haplotype identity, and calculate difference in number of polymorphic sites between haplotype groups

# cdsHaps.test <- copy(cdsHaps)
# cdsHaps.test[, haplotype := NULL]
# 
# # important to exclude fixed polymorphisms
# nonsyn_poly_sites <- cds_SNP_counts[type == 'N' & fixed_poly != 'fixed', pos]
# cdsHaps.test <- cdsHaps.test[pos %in% nonsyn_poly_sites]
# 
# diff.i <- vector()
# for(i in 1:10000){
#   # random ordering of haps
#   newHaps <- data.table(hapRank = 1:226, haplotype = sample(c(rep('H', 54), rep('h', 172)), size = 226, replace = F))
#   cdsHaps.test.i <- merge(cdsHaps.test, newHaps)
#   
#   # get allele counts by haplotype at each site
#   
#   cds_SNP_counts.i <- cdsHaps.test.i[, .(count = sum(allele)), by = .(pos, haplotype, allele0, allele1)]
#   cds_SNP_counts.i <- dcast(cds_SNP_counts.i, pos + allele0 + allele1 ~ haplotype, value.var = 'count')
#   
#   # make sure you double check final numbers when assigning these classes
#   cds_SNP_counts.i[H == 54 & h == 0, fixed_poly := 'fixed']
#   cds_SNP_counts.i[h > 0 & H < 54, fixed_poly := 'polymorphic_both']
#   cds_SNP_counts.i[h > 0 & H == 54, fixed_poly := 'polymorphic_h']
#   cds_SNP_counts.i[H < 54 & h == 0, fixed_poly := 'polymorphic_H']
#   
#   # calculate difference in permuted data set 
#   polyN.H.i <- cds_SNP_counts.i[fixed_poly == 'polymorphic_H', .N]
#   polyN.h.i <- cds_SNP_counts.i[fixed_poly == 'polymorphic_H', .N]
#   
#   diff.i[i] <- polyN.H.i - polyN.h.i
# }



# option to write file out, then samtools faidx the reference to get the cds sequence
# (had previously done, just read at top of script)
# fwrite(data.table(CDS[, paste0("NC_049911.1:", pos, "-", pos)]), 
#       file= "~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/TPPD_CDS_loc.txt",
#       quote = F, row.names = F, col.names = F)


