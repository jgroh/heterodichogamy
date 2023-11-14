library(ape)
library(data.table)
library(ggplot2)

gff <- read.gff("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/GCF_001411555.2_Walnut_2.0_genomic.gff")
setDT(gff)
TPP <- gff[grepl('XM_018957009.2', x = attributes, ignore.case = T) & type != 'mRNA', .(type, start, end)]

# read alignments

HJ1.0_TPPD <- rbindlist(lapply(list.files("/Users/Jeff/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments", 
           recursive = T, 
           pattern = '*HJ1-0_TPPD.txt',
           full.names = T), function(x){
             y <- fread(x)
             y[, qryGnom := gsub("/HJ1-0_TPPD.txt", "", gsub("/Users/Jeff/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_hJregCha_alignments/", "", x))]
             return(y)
           }
       ))
HJ1.0_TPPD[, alignmentID := NULL]

# make a qry seq entry for the reference, to reformat all sequences in the same column
HJ1.0_TPPD_hJregCha <- copy(HJ1.0_TPPD)
HJ1.0_TPPD_hJregCha <- HJ1.0_TPPD_hJregCha[qryGnom == 'HJregBNU']
HJ1.0_TPPD_hJregCha[, c("qryChr", "qryStart", "qryLength", "qrySeq", "qryGnom") := .(refChr, refStart, refLength, refSeq, 'hJregCha')]

HJ1.0_TPPD <- rbind(HJ1.0_TPPD, HJ1.0_TPPD_hJregCha)

# returns logical, TRUE if a site is variable 
filter_SNPs <- function(x) {
  counts <- table(x)
  
  if(length(counts) == 1){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
  
SNP_filter <- HJ1.0_TPPD[, .(filter_SNPs(qrySeq)), by = refPos]
SNP_positions <- SNP_filter[V1 == TRUE, refPos]

HJ1.0_TPPD_SNPs <- HJ1.0_TPPD[refPos %in% SNP_positions]

filter_singletons <- function(x){
  counts <- table(x)
  
  # filter biallelic 
  if(length(counts) == 2){
    if(min(counts) > 1){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
  
  # filter multi-allelic
  if(length(counts) > 2){
    if(length(unique(counts)) > 2){ # If there are only singletons present, there will be only two values of allele counts
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}

non_singletons <- HJ1.0_TPPD_SNPs[, .(non_singletons = filter_singletons(qrySeq)), by = refPos]
non_singleton_positions <- non_singletons[non_singletons == TRUE, refPos]

HJ1.0_TPPD_SNPs_fltd <- HJ1.0_TPPD_SNPs[refPos %in% non_singleton_positions]


HJ1.0_TPPD_SNPs_fltd[, qrySeq := factor(qrySeq, levels = c("-", "A", "T", 'C', "G"))]
HJ1.0_TPPD_SNPs_fltd[, rank := seq(1, .N), by = qryGnom]


HJ1.0_TPPD_SNPs_fltd[qryGnom == 'hJregCha', ID := 'h Jreg']
HJ1.0_TPPD_SNPs_fltd[qryGnom == 'HJregBNU', ID := 'H Jreg']
HJ1.0_TPPD_SNPs_fltd[qryGnom == 'hJcaliAlt', ID := 'h Jcal']
HJ1.0_TPPD_SNPs_fltd[qryGnom == 'HJcaliPrimary', ID := 'H Jcal']
HJ1.0_TPPD_SNPs_fltd[qryGnom == 'hJmanBNU', ID := 'h Jman']
HJ1.0_TPPD_SNPs_fltd[qryGnom == 'HJmanNFU', ID := 'H Jman']


ggplot(HJ1.0_TPPD_SNPs_fltd[refPos >= TPP[, min(start)] & refPos < TPP[, max(end)]], aes(x = rank, y = ID, fill = qrySeq, label = qrySeq)) + 
  geom_tile(width = 1, height = 1) + 
  #geom_text(size = 3, color = 'black',show.legend = FALSE, aes(label = ifelse(base != "-", as.character(base), ""))) +
  geom_text(size = 3, color = 'black',show.legend = FALSE, aes(as.character(qrySeq))) +
  theme_classic() + 
  scale_fill_manual(values = c("gray", 'blue', 'yellow', 'green', 'red')) +
  labs(x = '', y = '') + 
  #scale_x_discrete(labels = sort(unique(HJ1.0_TPPD_SNPs_fltd$refPos)), position = 'top') + 
  theme(aspect.ratio = 0.1,
        legend.position = 'none', 
        plot.margin = unit(c(0,1,0,0), 'cm'),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 40, size = 7, hjust = 0),
        axis.ticks.length.x = unit(0, 'cm'))  






n_unique <-  HJ1.0_TPPD[, .(unique_values = length(unique(qrySeq))), by = refPos]
non_singleton_positions <- n_unique[unique_values >= 2, refPos]

HJ1.0_TPPD

map <- fread("~/workspace/heterodichogamy/H_locus_structure/TPPD/gene_w_introns/out.map")
ped <- fread("~/workspace/heterodichogamy/H_locus_structure/TPPD/gene_w_introns/out.ped")
pos <- map[, .(pos = V4)]

nuc <- ped[,-c(2:6)]
setnames(nuc, "V1", "ID")

# each column is duplicated immediately after because it thinks they're diploid
selectcols <- c("ID", paste0("V", seq(7, 110, by = 2)))
nuc <- nuc[, ..selectcols]

setnames(nuc, paste0("V", seq(7, 110, by = 2)), as.character(1:52))

# to long, clean up for plotting
nuc <- melt(nuc, id.vars = "ID", variable.name = "rank", value.name = "base")
nuc[, pos := pos[, pos], by = ID]

nuc[ID == 'hJregCha.bam', ID := 'h Jreg']
nuc[ID == 'HJregBNU.bam', ID := 'H Jreg']
nuc[ID == 'hJcaliAlt.bam', ID := 'h Jcal']
nuc[ID == 'HJcaliPrimary.bam', ID := 'H Jcal']
nuc[ID == 'hJmanBNU.bam', ID := 'h Jman']
nuc[ID == 'HJmanNFU.bam', ID := 'H Jman']


nuc[, base := factor(base, levels = c("0", "A", "T", 'C', "G"))]
ggplot(nuc, aes(x = rank, y = ID, fill = base, label = base)) + 
  geom_tile(width = 1, height = 1) + 
  geom_text(size = 3, color = 'black',show.legend = FALSE, aes(label = ifelse(base != "0", as.character(base), ""))) +
  theme_classic() + 
  scale_fill_manual(values = c("gray", 'blue', 'yellow', 'green', 'red')) +
  labs(x = '', y = '') + 
  scale_x_discrete(labels = sort(unique(nuc$pos)), position = 'top') + 
  theme(aspect.ratio = 0.1,
        legend.position = 'none', 
        plot.margin = unit(c(0,1,0,0), 'cm'),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 40, size = 7, hjust = 0),
        axis.ticks.length.x = unit(0, 'cm'))  


plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
#lines(x = c(0, 1), y = c(0, 0))
lines(x = c(0, 1), y = c(.4, .4))

pos_scld <- unlist(pos) 
pos_scld <- (pos_scld - min(TPP[, start]))/(max(TPP[, end]) - min(TPP[, start]))
for(i in 1:52){
  lines(x = c(i/52, pos_scld[i]), y= c(0.5, 0.4), lwd = 0.5)
}

TPP[, start_scld := (start - min(start))/(max(end) - min(start))]
TPP[, end_scld := (end - min(start))/(max(end) - min(start))]

for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'exon']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.35, 0.4, 0.4, 0.35), col = 'gray')
  }
}
for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'CDS']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.35, 0.4, 0.4, 0.35), col = 'salmon')
  }
}




