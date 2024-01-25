library(data.table)

# This script reformats alignment result between pecan 'Pawnee' and J. regia V2 assemblies

maf <- '~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments/CillPaw/Carya_chr4.maf'

# read alignment. Use sep value that's not in file to force fread to behave like readLines
mafLines <- fread(maf, sep = "?", skip = 1,
                  header = F)


# ===== format maf alignment into data.table =====

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

# *swqp reference and query* so that pecan is the reference here
x <- c(2,1)
for(i in 1:((nrow(mafDT) - 1)/2)){ x <- c(x, c(2,1) + 2*i) }
mafDT <- mafDT[x]

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


# ====== Tally Substitutions ======
# assign logical to substitutions
seqDT[, sub := ifelse(refSeq != qrySeq, T, F)]

# gaps in the query sequence are treated as NA 
# (we account for the total number of non-gapped bps in the calc of Dxy)
seqDT[qrySeq == "-", sub := NA]

# check total divergence 
# seqDT[, sum(sub, na.rm=T)/length(which(!is.na(sub)))]




