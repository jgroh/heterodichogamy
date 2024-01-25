pos <- fread("~/workspace/heterodichogamy/pecan/WGS_data/homozygotes_chr4_6.4-7.05Mb.012.pos")
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/homozygotes_chr4_6.4-7.05Mb.012")
indv1 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/homozygotes_chr4_6.4-7.05Mb.012.indv", header = F, col.names = c("ID"))
ids <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_geno.txt")
indv2 <- merge(ids, indv)
indv <- indv2[indv1[,ID]]

setnames(geno, c("ID", pos[, V2]))

HH <- geno[which(indv$genotype == 'HH')]
hh <- geno[which(indv$genotype == 'hh')]

index_hh <- apply(hh[,-1], MARGIN = 2, function(x){all(x %in% c(2, -1))})
index_HH <- apply(HH[,-1], MARGIN = 2, function(x){all(x %in% c(0, -1))})

pos <- intersect(names(index_hh[index_hh]), names(index_HH[index_HH]))

output <- data.table(chr = 'CM031828.1', pos = pos)

fwrite(output, file = '~/workspace/heterodichogamy/pecan/WGS_data/Hloc_fixed_CDS_SNPs.txt', quote = F, col.names = F, row.names = F, sep = "\t")
