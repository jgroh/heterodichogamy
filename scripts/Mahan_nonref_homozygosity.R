library(data.table)
library(ggplot2)

pos <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Mahan_CDS.012.pos")
pos <- pos[, V2]
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Mahan_CDS.012")
geno <- melt(geno, value.name = 'genotype')
geno <- geno[-1]
geno[, variable := NULL]

geno[, pos := pos]

# looks as expected
ggplot(geno[genotype != -1 & pos > 6.42e6 & pos < 6.7e6], aes(x = pos, y = genotype)) + 
  geom_point() + geom_smooth(method = 'loess', span = 0.5) + 
  labs(x = "Pawnee chr 4", y = "Number of non-reference alleles in Mahan") +
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

# select sites
output <- geno[genotype == 2 & pos > 6461404 & pos < 6658636]
output[, genotype := NULL]
output[, chr := 'CM031812.1']

fwrite(x = output[, .(chr, pos)], 
       file = '~/workspace/heterodichogamy/pecan/WGS_data/Mahan_homozygous_nonref_sites.txt',
       col.names = F, row.names = F, quote = F, sep = '\t')
