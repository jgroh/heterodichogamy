library(data.table)

vcfcols <- fread("~/workspace/heterodichogamy/regia/calls/vcfcols.txt", col.names = 'ID', header = F)
vcfcols[10:nrow(vcfcols)]

vcfcols1 <- vcfcols[-c(1:9)]

# merge with genotypes
geno <- fread("~/workspace/heterodichogamy/regia/sra_samples_geno.txt")
  
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("ID", "phenotype", 'species'), header = F)
pheno[species == 'regia' & phenotype == 'protandrous', genotype := 'gg']
pheno[species == 'regia' & phenotype == 'protogynous' & ID != 'JG0026', genotype := 'Gg']
pheno[species == 'regia' & phenotype == 'protogynous' & ID == 'JG0026', genotype := 'GG']


regia_genotypes <- rbind(geno[, .(ID=run, genotype)], pheno[species == 'regia', .(ID, genotype)])


d1 <- merge(vcfcols1, regia_genotypes)


vcfcols1[, ID]

vcfcols2 <- d1[vcfcols1[, ID]]

# change genotypes and position depending on reference used
writegt <- function(gt){
  if(gt == 'gg'){
    x <- "0/0:0,135,255:45:45,0:1,6.31546e-15,3.15319e-28:99"
  } else if(gt == 'Gg'){
    x <- "0/1:255,0,255:32:19,13:3.12468e-26,1,8.00083e-27:99"
  } else if(gt == 'GG'){
    x <- "1/1:255,75,0:25:0,25:1.37819e-25,1.32033e-07,1:68"
  }
}

vcfcols2[, vcfgeno := writegt(genotype), by = ID]

part1 <- "NC_049911.1\t31884475\t.\tG\tA\t2827.3\t.\t.\tGT:PL:DP:AD:GP:GQ\t"

part2 <- paste0(vcfcols2[, unlist(vcfgeno)], collapse = '\t')

write(paste0(part1, part2), file = '~/workspace/heterodichogamy/regia/calls/fake_SNP_line.txt')
