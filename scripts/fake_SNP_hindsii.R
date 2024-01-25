library(data.table)

vcfcols <- fread("~/workspace/heterodichogamy/hindsii/vcfcols.txt", col.names = 'ID', header = F)
vcfcols[10:nrow(vcfcols)]

vcfcols1 <- vcfcols[-c(1:9)]

# merge with genotypes
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", header = F, col.names = c("ID", "phenotype", "species"))

pheno[species == 'hindsii' & phenotype == 'protogynous', genotype := 'Gg']
pheno[species == 'hindsii' & phenotype == 'protandrous', genotype := 'gg']

d1 <- merge(vcfcols1, pheno)

d1[ID == 'JHIN_PC_102', genotype := 'Gg']
d1[ID == 'JHIN_PC_052', genotype := 'gg']


vcfcols1[, ID]

vcfcols2 <- d1[vcfcols1[, ID]]


writegt <- function(gt){
  if(gt == 'Gg'){
    x <- "0/1:0,117,255:39:39,0:1,3.84607e-13,2.93747e-28:99"
  } else if(gt == 'gg'){
    x <- "1/1:255,111,0:37:0,37:5.09415e-27,6.37626e-12,1:99"
  }
}

vcfcols2[, vcfgeno := writegt(genotype), by = ID]

part1 <- "JAKSXK010000007.1\t31370000\t.\tG\tA\t2827.3\t.\t.\tGT:PL:DP:AD:GP:GQ\t"

part2 <- paste0(vcfcols2[, unlist(vcfgeno)], collapse = '\t')

write(paste0(part1, part2), file = '~/workspace/heterodichogamy/hindsii/fake_SNP_line.txt')
