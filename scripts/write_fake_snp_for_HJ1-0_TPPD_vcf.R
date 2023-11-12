library(data.table)

vcfcols <- fread("~/workspace/heterodichogamy/regia/calls/vcfcols.txt", col.names = 'run', header = F)
vcfcols[10:nrow(vcfcols)]

vcfcols1 <- vcfcols[-c(1:9)]

# merge with genotypes
sra_geno <- fread("~/workspace/heterodichogamy/regia/sra_samples_geno.txt")
sra_geno[, c('ID', 'loc') := NULL]

pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", header = F, col.names = c("ID", "phenotype", "species"))
pheno[ID != 'JG0026' & phenotype == 'protogynous' & species == 'regia', genotype := 'Hh']
pheno[ID == 'JG0026', genotype := 'HH']
pheno[phenotype=='protandrous', genotype := 'hh']
reg_geno <- pheno[species == 'regia', .(run = ID, genotype)]

genotypes <- rbind(reg_geno, sra_geno)

d1 <- merge(vcfcols1, genotypes)

vcfcols1[, run]

vcfcols2 <- d1[vcfcols1[, run]]


writegt <- function(gt){
  if(gt == 'HH'){
    x <- "1/1:0,117,255:39:39,0:1,3.84607e-13,2.93747e-28:99"
  } else if(gt == 'Hh'){
    x <- "0/1:255,0,255:46:28,18:9.02354e-27,1,2.77053e-26:99"
  } else if(gt == 'hh'){
    x <- "0/0:255,111,0:37:0,37:5.09415e-27,6.37626e-12,1:99"
  }
}

vcfcols2[, vcfgeno := writegt(genotype), by = run]

part1 <- "NC_049911.1\t31875000\t.\tG\tA\t2827.3\t.\t.\tGT:PL:DP:AD:GP:GQ\t"

part2 <- paste0(vcfcols2[, unlist(vcfgeno)], collapse = '\t')

write(paste0(part1, part2), file = '~/workspace/heterodichogamy/regia/calls/fake_snp_hJ2.txt')
