library(data.table)
library(readxl)

fam.1 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars_old.fam")

pheno_my <- setDT(read_excel("~/workspace/heterodichogamy/Manuscript/Supplementary_tables.xlsx", sheet = 'NovaSeq'))
pheno_sra <- fread("~/workspace/heterodichogamy/pecan/WGS_data/sra_samples_subset.txt")
pheno <- rbind(pheno_sra[, .(V2 = run, V6 = phenotype)],pheno_my[, .(V2 = ID, V6 = phenotype)])

fam.1[, V6 := NULL]

fam <- merge(fam.1, pheno)

fam[V6 == 'protandrous', V6 := 0]
fam[V6 == 'protogynous', V6 := 1]


fwrite(fam[fam.1[,V2]], file = '~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.fam', sep = "\t", col.names = F, row.names = F, quote = F)
