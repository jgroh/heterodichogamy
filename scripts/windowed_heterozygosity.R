library(data.table)
library(ggplot2)
library(magrittr)
library(ggpubr)


# ----- read pecan data 
het <- fread("~/workspace/heterodichogamy/pecan/WGS_data/pixy_pi.txt")
smpl <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_reassigned.tsv")
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c("run", 'genotype'))

het <- merge(het, smpl[, .(run, pop= type, phenotype)])
het <- merge(het, geno, by = 'run')
het[, pos := (window_pos_1 + window_pos_2)/2]

# ----- read Carya spp data
het_spp <- fread("~/workspace/heterodichogamy/Carya_spp/pixy_pi.txt")
het_spp[, pos := (window_pos_1 + window_pos_2)/2]

ENA <- c("aquatica", "cordiformis", "floridana", "glara", "laciniosa", 
         "myristiciformis", "ovata", "palmeri", "texana", "tomentosa")
EA <- c("cathayensis", "dabieshanensis", "hunanensis", 'kweichowensis', 'poilanei', "tonkinensis")

het_spp[pop %in% ENA, clade := 'ENA']
het_spp[pop %in% EA, clade := 'EA']
# ----- plot pecan heterozygosity

pecanplot <- ggplot(het[window_pos_1 > 6.35e6 & window_pos_2 < 6.75e6], 
       aes(x = pos, y = avg_pi, group = pop, color = genotype)) + 
  geom_point()  + 
  geom_smooth(method = 'loess', se = F, span = 0.5) + 
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  scale_x_continuous(breaks = seq(6.35e6, 6.75e6, length.out = 6), labels = seq(6.35, 6.75, length.out = 6))  + # if bottom
  theme_classic() + 
  labs(x = 'Pawnee chr 4 (Mb)', 
       y = 'Heterozygosity per exon', 
       color = '',
       title = 'Pecan') +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),

        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13), 
        #legend.position = c(0.1, 0.8)
        )
  

# ----- plot Carya spp heterozygosity

ENAplot <- ggplot(het_spp[clade == 'ENA' & window_pos_1 > 6.35e6 & window_pos_2 < 6.75e6], 
       aes(x = pos, y = avg_pi, group = pop, color = pop)) + 
  #facet_wrap(~clade) +
  scale_color_brewer(type = 'qual', palette = 'Paired') +
  #Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  geom_point()  + 
  
  geom_smooth(method = 'loess', se = F, span = 0.5) + 
  #scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  scale_x_continuous(breaks = seq(6.35e6, 6.75e6, length.out = 6), labels = seq(6.35, 6.75, length.out = 6))  + # if bottom
  theme_classic() + 
  labs(x = 'Pawnee chr 4 (Mb)', 
       y = 'Heterozygosity per exon', color = '', title = 'North America clade') +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))

#ENAplot


EAplot <- ggplot(het_spp[clade == 'EA' & window_pos_1 > 6.35e6 & window_pos_2 < 6.75e6], 
       aes(x = pos, y = avg_pi, group = pop, color = pop)) + 
  #facet_wrap(~clade) +
  scale_color_brewer(type = 'qual', palette = 'Paired') +
  #Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  geom_point()  + 
  
  geom_smooth(method = 'loess', se = F, span = 0.5) + 
  #scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  scale_x_continuous(breaks = seq(6.35e6, 6.75e6, length.out = 6), labels = seq(6.35, 6.75, length.out = 6))  + # if bottom
  theme_classic() + 
  labs(x = 'Pawnee chr 4 (Mb)', 
       y = 'Heterozygosity per exon', color = '', title = 'East Asian clade') +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))
  


plot_grid(pecanplot, ENAplot, EAplot, align = 'v')





# ----- first check expected pattern of heterozygosity for known genotypes in pecan -----
ids <- fread("~/workspace/heterodichogamy/pecan/WGS_data/out.012.indv", header = F)
pos <- fread("~/workspace/heterodichogamy/pecan/WGS_data/out.012.pos", header = F)
het0 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/out.012", header = F)
bed <- fread("~/workspace/heterodichogamy/Carya_spp/Pawnee_chr4_genes.bed", col.names = c("chr", 'start', 'end', 'gene', 'idk', 'strand'))
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c("INDV", "geno"))

# reformat
het0[, V1 := ids$V1]
setnames(het0, c("ID", pos$V2))
het1 <- melt(het0, id.vars = 'ID', variable.name = "pos", value.name = 'genotype')
het1[, pos := as.numeric(as.character(pos))]
head(het1)

het1

het1[, window := cut(pos, breaks = seq(min(pos), max(pos)+5000, by = 5000), labels = seq(min(pos), max(pos), by = 5000), include.lowest =T), by = ID]
het1[, het := ifelse(genotype == 1, 1, 0)]
het1[genotype == -1, het := NA]
het2 <- het1[, .(het1kb = mean(het, na.rm = T)), by = .(ID, window)]
het2[, window := as.numeric(as.character((window)))]

het <- merge(het2, geno[, .(ID = INDV, geno)])

het[window > 6.42e6 & window < 6.7e6, mean(het1kb), by = .(ID, geno)] %>%
  ggplot(., aes(x = geno, y = V1)) + geom_point()

ggplot(het[window > 6.42e6 & window < 6.7e6], aes(x = window, y = het1kb, color = geno)) + 
  geom_point() + 
  scale_color_manual(values = c("lightgray", "salmon", "black")) + 
  theme_classic() + 
  geom_smooth(aes(group = ID), method = 'loess', span = 0.3, se = F)

#SRR15911540 looks like a heterozygote despite not having structural variant - weird. 

# get list of positions that fall within genes
gene_pos <- NULL
for(p in unique(het1$pos)){
  if(any(bed[start > 5e6 & end < 8e6, p >= start & p <= end])){
    gene_pos <- c(gene_pos, p)
  }
}







# ------------ Carya spp -----

ids.1 <- fread("~/workspace/heterodichogamy/Carya_spp/out.012.indv", header = F)
pos.1 <- fread("~/workspace/heterodichogamy/Carya_spp/out.012.pos", header = F)
het0.1 <- fread("~/workspace/heterodichogamy/Carya_spp/out.012", header = F)


# reformat
het0.1[, V1 := ids.1$V1]
setnames(het0.1, c("ID", pos.1$V2))
het1.1 <- melt(het0.1, id.vars = 'ID', variable.name = "pos", value.name = 'genotype')
het1.1[, pos := as.numeric(as.character(pos))]
head(het1.1)

het1.1[, window := cut(pos, breaks = seq(min(pos), max(pos)+20000, by = 20000), labels = seq(min(pos), max(pos), by = 20000), include.lowest =T), by = ID]
het1.1[, het := ifelse(genotype == 1, 1, 0)]
het1.1[genotype == -1, het := NA]
het2.1 <- het1.1[, .(het1kb = mean(het, na.rm = T)), by = .(ID, window)]
het2.1[, window := as.numeric(as.character((window)))]


ggplot(het2.1[window > 6.3e6 & window < 7e6], aes(x = window, y = het1kb, color = ID)) + 
  geom_point() + 
  theme_classic() +  geom_line()
 # geom_smooth(aes(group = ID), method = 'loess', span = 0.2, se = F)

possible_hets <- c('SRR6804841', 'SRR6804842','SRR6804843', 'SRR6804844','SRR6804845', 'SRR6804847')
ggplot(het2.1[window > 6.3e6 & window < 6.8e6 & ID %in% possible_hets], aes(x = window, y = het1kb, color = ID)) + 
  geom_point() + 
  theme_classic() + 
  geom_smooth(aes(group = ID), method = 'loess', span = 0.4, se = F)





