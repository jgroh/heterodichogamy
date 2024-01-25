library(data.table)

# Calculate Dxy from multiple alignment format output from Anchorwave. 
# ! This script is for a specific alignment and would need to be edited for general use. 
mafLines <- fread("~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_Csin/Ccat/alignment.maf", sep = "?", header =F, skip = 1)

chrom <- "Chr04"


# ------- read coding regions -----
cds <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Csin_cds_coords.txt", col.names = c("chr", "start", "end"))
cds <- cds[chr == 'Chr04']
setkey(cds, start)

# ===== format maf alignment into data.table =====

# subset to alignment blocks
mafLines <- mafLines[grep("^s", V1)]

# set column names
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5
mafDT <- mafLines[, tstrsplit(V1, "\\s+")]
setnames(mafDT, c("s", "chr", "start0", "length", "strand", "srcLength", "seq"))

# subset to focal chromosome. Find which rows have the chromosome ID, and also take the next row which has the aligned chromosome
# in this case the two aligned genomes have same chromosome ID, so how to select index

indices <- which(mafDT$chr == chrom) # indices for ref chrom
indices <- indices[indices %% 2 == 1]
indices <- c(indices, indices + 1) # indices for qry chrom
indices <- indices[order(indices)]
mafDT <- mafDT[indices]
    
# add grouping variable for alignment id
mafDT[, alignmentID := rep(1:(.N / 2), each = 2)]
# mafDT[, c(1:6,8)]

# subset to focal columns
mafDT <- mafDT[, .(alignmentID, chr, start0 = as.numeric(start0), length = as.numeric(length), seq)]
# mafDT[, c(1:4)]
mafDT0 <- copy(mafDT)

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
  
  # For calculating Dxy against a reference, 
  # sites in the alignment that are absent in the reference are not informative, so we disregard these. 
  outDT <- outDT[refSeq != "-"]
  
  # set coordinates of reference. 
  outDT[, refPos := seq_len(.N)]
  outDT[, refPos := refPos + refStart - 1]
  outDT
}


seqDT <- mafDT[, tr(DT = .SD, by = .BY), by = alignmentID]



# ====== Tally Substitutions ======

# assign logical to substitutions
seqDT[, sub := ifelse(refSeq != qrySeq, T, F)]

# gaps in the query sequence are treated as NA 
# (we account for the total number of non-gapped bps in the calc of Dxy)
seqDT[qrySeq == "-", sub := NA]


# ===== Assign each bp to its containing window =====
round_up <- function(x) {
  result <- ceiling(x / window_size) * window_size
  return(result)
}

round_down <- function(x) {
  result <- floor(x / window_size) * window_size
  return(result)
}

window_size <- 1e3
seqDT[, window := cut(refPos, 
                      breaks = seq(from = round_down(refStart[1]), to = round_up(refStart[1] + refLength[1]), by = window_size) , 
                      include.lowest=T, 
                      labels = seq(from = round_down(refStart[1]), 
                                   to = round_up(refStart[1] + refLength[1]) - window_size, by = window_size) + 0.5*window_size), 
      by = alignmentID]


# ===== compute Dxy in windows =====
Dxy <- seqDT[, .(Dxy_num = sum(sub, na.rm=T), 
                 Dxy_denom = sum(!is.na(sub))), by = .(alignmentID, refChr, qryChr, window)]
Dxy[, Dxy := Dxy_num/Dxy_denom]
Dxy[Dxy_num == 0 & Dxy_denom == 0, Dxy := NA]

#Dxy[, refGnom := ref]
#Dxy[, qryGnom := qry]

Dxy[, window := as.numeric(as.character(window))]

ggplot(Dxy[window > 5600000 & window < 6900000], aes(x = window, y = Dxy)) + 
  geom_point() + 
  geom_smooth(method = 'loess') + 
  labs(y = '', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. cathayensis') + # if bottom
  scale_x_continuous(breaks = seq(5.7e6, 6.7e6, length.out = 6), labels = seq(5.7, 6.7, length.out = 6) )  +
  theme_classic() + 
  theme(aspect.ratio = .2,
        #plot.margin = margin(0,0,0,0, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        #legend.position = c(0.15, 0.85),
        plot.title = element_text(face = 'italic'),
        legend.position = 'none',
        legend.key.size = unit(0.3, 'cm'),
  )




# ===== compute Dxy for coding regions ===== (code not working, used different approach)
# all_sub <- seqDT[, .(refPos, sub)]
# 
# sub_sub <- all_sub[refPos > 5000000 & refPos < 8000000]
# cds_sub <- cds[start > 5000000 & end < 8000000]


# all_sub_fltd <- all_sub[cds, on = .(refPos >= start, refPos <= end), nomatch = 0L]
# #
# result <- cds[all_sub, 
#               .(dxy_num = sum(sub, na.rm = TRUE), 
#                 dxy_denom = sum(!is.na(sub))), 
#               on = .(chr = refChr, start <= refPos, end >= refPos), 
#               by = .EACHI]
# 
# 
# 
# cds_dxy <- cds[seqDT, .(dxy_num = sum(sub, na.rm = TRUE), 
#                  dxy_denom = sum(!is.na(sub))), 
#             on = .(start <= refPos, end >= refPos), 
#             by = .EACHI]


# seqDT_CDS <- seqDT[cds, on = .(refPos >= start, refPos <= end), nomatch = 0L]


# for(i in 1:nrow(cds_sub)){
#   dxy_num <- sub_sub[refPos >= cds_sub[i, start] & refPos <= cds_sub[i, end], sum(sub, na.rm = T)]
#   dxy_denom <- sub_sub[refPos >= cds_sub[i, start] & refPos <= cds_sub[i, end], sum(!is.na(sub))]
#   cds_sub[i, dxy := dxy_num/dxy_denom]
# }
#   
# 
# 
# 
# ggplot(cds_sub[start > 5.6e6 & end < 6.9e6], aes(x = (start + end)/2, y = dxy) ) + geom_point()






# ===== Output =====
#fwrite(Dxy, file = "", quote = F, sep = "\t", row.names = F, na = "NA")
