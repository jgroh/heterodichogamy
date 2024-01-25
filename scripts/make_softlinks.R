ln <- fread("~/workspace/heterodichogamy/newdata/softlinks2.txt", header = F)

ln[,unique(substr(V2, 1, 4))]

ln[substr(V2, 1, 4) == 'CILL', V3 := paste0('/home/jgroh/heterodichogamy/pecan/WGS_data/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'COVA', V3 := paste0('/home/jgroh/heterodichogamy/Carya_spp/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JAIL', V3 := paste0('/home/jgroh/heterodichogamy/ailantifolia/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JCAL', V3 := paste0('/home/jgroh/heterodichogamy/californica/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JCAT', V3 := paste0('/home/jgroh/heterodichogamy/cathayensis/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JCIN', V3 := paste0('/home/jgroh/heterodichogamy/cinerea/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JHIN', V3 := paste0('/home/jgroh/heterodichogamy/hindsii/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'JMAJ', V3 := paste0('/home/jgroh/heterodichogamy/major/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'PMAC', V3 := paste0('/home/jgroh/heterodichogamy/Pterocarya_macroptera/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'PRHO', V3 := paste0('/home/jgroh/heterodichogamy/Pterocarya_rhoifolia/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'PSTE', V3 := paste0('/home/jgroh/heterodichogamy/Pterocarya_stenoptera/fastq/', V2, ".fq.gz")]
ln[substr(V2, 1, 4) == 'PSTR', V3 := paste0('/home/jgroh/heterodichogamy/Platycarya_strobilaceae/fastq/', V2, ".fq.gz")]

ln[, V2 := NULL]

ln[, V1 := paste0("/home/jgroh/heterodichogamy/newdata/fastq/jzb7z0gzgs/Un_DTSA841/Project_CLDV_Nova974P_Vik/", V1)]
fwrite(ln, file = '~/workspace/heterodichogamy/newdata/softlinks.txt', col.names = F, row.names = F, quote = F, sep = "\t")
