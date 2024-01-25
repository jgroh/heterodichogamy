library(data.table)

vcfcols <- fread("~/workspace/heterodichogamy/pecan/WGS_data/vcfcols.txt", col.names = 'run', header = F)
vcfcols[10:nrow(vcfcols)]

vcfcols1 <- vcfcols[-c(1:9)]

# merge with genotypes
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")
geno <- geno[, .(run = ID, genotype)]


d1 <- merge(vcfcols1, geno)

vcfcols1[, run]

vcfcols2 <- d1[vcfcols1[, run]]


writegt <- function(gt){
  if(gt == 'HH'){
    x <- "0/0:0,117,255:39:39,0:1,3.84607e-13,2.93747e-28:99"
  } else if(gt == 'Hh'){
    x <- "0/1:255,0,255:46:28,18:9.02354e-27,1,2.77053e-26:99"
  } else if(gt == 'hh'){
    x <- "1/1:255,111,0:37:0,37:5.09415e-27,6.37626e-12,1:99"
  }
}

vcfcols2[, vcfgeno := writegt(genotype), by = run]

part1 <- "CM031828.1\t6700000\t.\tG\tA\t2827.3\t.\t.\tGT:PL:DP:AD:GP:GQ\t"

part2 <- paste0(vcfcols2[, unlist(vcfgeno)], collapse = '\t')

write(paste0(part1, part2), file = '~/workspace/heterodichogamy/pecan/WGS_data/fake_SNP_line.txt')
