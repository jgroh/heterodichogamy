library(data.table)
library(readxl)

fam1 <- fread("~/workspace/heterodichogamy/hindsii/calls/Jcali_alt/Putah2JcaliAlt_old.fam")
fam1

pheno <- read_excel("~/workspace/heterodichogamy/Manuscript/Supplementary_tables.xlsx", sheet = "NovaSeq")
setDT(pheno)


fam2 <- merge(fam1[, .(V1, V2, V3, V4, V5)], pheno[, .(V2 = ID, V6 = phenotype)])

fam2 <- fam2[, .(V1, V2, V3, V4, V5, V6)]
      
fam2[V6 == 'protandrous', V6 := 0]
fam2[V6 == 'protogynous', V6 := 1]

fam2[V6 %in% 0:1]

fwrite(fam2[V6 %in% 0:1], file = '~/workspace/heterodichogamy/hindsii/calls/Jcali_alt/Putah2JcaliAlt.fam',
       col.names = F, row.names = F, quote = F, sep = "\t")
