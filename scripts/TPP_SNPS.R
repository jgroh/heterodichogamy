library(ape)
library(data.table)
library(ggplot2)

gff <- read.gff("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/GCF_001411555.2_Walnut_2.0_genomic.gff")
setDT(gff)
TPP <- gff[grepl('XM_018957009.2', x = attributes, ignore.case = T) & type != 'mRNA', .(type, start, end)]

map <- fread("~/workspace/heterodichogamy/H_locus_structure/HJ1-0_TPPD/out.map")
ped <- fread("~/workspace/heterodichogamy/H_locus_structure/HJ1-0_TPPD/out.ped")
pos <- map[, .(pos = V4)]

nuc <- ped[,-c(2:6)]
setnames(nuc, "V1", "ID")

# each column is duplicated immediately after because it thinks they're diploid
# how many SNPs?
mx <- max( as.numeric( gsub("V", "", names(nuc)[-1] ) ) )
SNP_nx2 <- mx - 7 + 1
nSNP <- SNP_nx2/2

selectcols <- c("ID", paste0("V", seq(7, mx, by = 2)))
nuc <- nuc[, ..selectcols]

setnames(nuc, paste0("V", seq(7, mx, by = 2)), as.character(1:(SNP_nx2/2)))

# to long, clean up for plotting
nuc <- melt(nuc, id.vars = "ID", variable.name = "rank", value.name = "base")
nuc[, pos := pos[, pos], by = ID]

nuc[ID == 'PsteBNU.bam', ID := 'Pste']
nuc[ID == 'hJregCha.bam', ID := 'h Jreg']
nuc[ID == 'HJregBNU.bam', ID := 'H Jreg']
nuc[ID == 'hJcalAlt.bam', ID := 'h Jcal']
nuc[ID == 'HJcalPrimary.bam', ID := 'H Jcal']
nuc[ID == 'hJmanBNU.bam', ID := 'h Jman']
nuc[ID == 'HJmanNFU.bam', ID := 'H Jman']
nuc[ID == 'hJnigBNU.bam', ID := 'h Jnig']
nuc[ID == 'HJmic.bam', ID := 'H Jmic']


nuc[, base := factor(base, levels = c("0", "A", "T", 'C', "G"))]
nuc[, ID := factor(ID, levels = rev(c("Pste", "H Jreg", "h Jreg", "H Jman", "h Jman", 'H Jcal', 'h Jcal', 'H Jmic', 'h Jnig')))]

ggplot(nuc, aes(x = rank, y = ID, fill = base, label = base)) + 
  geom_tile(width = 1, height = 1) + 
  geom_text(size = 3, color = 'black',show.legend = FALSE, aes(label = ifelse(base != "0", as.character(base), ""))) +
  theme_classic() + 
  scale_fill_manual(values = c("gray", 'blue', 'yellow', 'green', 'red')) +
  labs(x = '', y = '') + 
  scale_x_discrete(labels = sort(unique(nuc$pos)), position = 'top') + 
  theme(aspect.ratio = 0.1,
        legend.position = 'none', 
        axis.line.x = element_line(size = 0),
        plot.margin = unit(c(0,1,0,0), 'cm'),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 40, size = 6, hjust = 0),
        axis.ticks.length.x = unit(0, 'cm'))  


plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
#lines(x = c(0, 1), y = c(0, 0))
lines(x = c(0, 1), y = c(.4, .4))

pos_scld <- unlist(pos) 
pos_scld <- (pos_scld - min(pos))/(max(TPP[, end]) - min(pos))

for(i in 1:nSNP){
  lines(x = c(i/nSNP, pos_scld[i]), y= c(0.5, 0.4), lwd = 0.5)
}

TPP[, start_scld := (start - min(pos))/(max(end) - min(pos))]
TPP[, end_scld := (end - min(pos))/(max(end) - min(pos))]

for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'exon']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.35, 0.4, 0.4, 0.35), col = 'gray')
  }
}
for(i in 1:nrow(TPP)){
  if(TPP[i, type == 'CDS']){
    polygon(x= c(TPP[i, start_scld], TPP[i, start_scld], TPP[i, end_scld], TPP[i, end_scld]),
            y = c(0.35, 0.4, 0.4, 0.35), col = 'goldenrod3')
  }
}




