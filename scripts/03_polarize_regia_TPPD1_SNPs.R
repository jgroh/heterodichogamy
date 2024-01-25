library(data.table)


# polarize fixed SNPs in TPPD1 between J. regia G-locus haplotypes


maf <- '~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments/CillPaw/chr11.maf'

# read alignment. Use sep value that's not in file to force fread to behave like readLines
mafLines <- fread(maf, sep = "?", skip = 1,
                  header = F)

# subset to alignment blocks
mafLines <- mafLines[grep("^s", V1)]

# set column names
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5
mafDT <- mafLines[, tstrsplit(V1, "\\s+")]
setnames(mafDT, c("s", "chr", "start0", "length", "strand", "srcLength", "seq"))


# add grouping variable for alignment id
mafDT[, alignmentID := rep(1:(.N / 2), each = 2)]
# mafDT[, c(1:6,8)]

# subset to focal columns
mafDT <- mafDT[, .(alignmentID, chr, start0 = as.numeric(start0), length = as.numeric(length), seq)]
# mafDT[, c(1:4)]


# function that takes alignment block and returns data table with aligned seqs in long format. First row in alignment block is ref, 2nd is qry.
tr <- function(DT, by){
  outDT <- data.table(alignmentID = by,
                      refChr = DT[1, chr], 
                      qryChr = DT[2, chr], 
                      refStart = DT[1, start0], 
                      qryStart = DT[2, start0],
                      refLength = DT[1, length], 
                      qryLength = DT[2, length],
                      refSeq = unlist(strsplit(DT[1, seq], "")),
                      qrySeq = unlist(strsplit(DT[2, seq], ""))
  )
  
  # ignore sites in the alignment that are absent in the reference so that coordinates can be earily set for reference in next step
  outDT <- outDT[refSeq != "-"]
  
  # set coordinates of reference. 
  outDT[, refPos := seq_len(.N)]
  outDT[, refPos := refPos + refStart]
  outDT
}

seqDT <- mafDT[, tr(DT = .SD, by = .BY), by = alignmentID]

seqDT

TPPD1 <- seqDT[refPos >= 31884478 & refPos <= 31887006]


TPPD1



