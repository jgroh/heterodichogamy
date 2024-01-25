library(data.table)
library(ggplot2)
library(ggpubr)

# ---- read pecan -----
d <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Mahan_homozygous_nonref_SNPs_pi.txt")

smpl <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_reassigned.tsv")
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c("run", 'genotype'))

d <- merge(d, smpl[, .(run, pop= type, phenotype)])
d <- merge(d, geno, by = 'run')

# calculate heterozygosity at SNPs
het <- d[, .(p = mean(avg_pi, na.rm=T)), by = run]

# calculate standard error
N <- d[no_sites == 1, .N, by = run]

het <- merge(N, het)
het <- merge(het, geno)

# exclude Mahan
setkey(het, genotype)
lv <- rev(unlist(het[, run]))
het[, run := factor(run, levels = lv)]

# -----pecan, all samples ----- 
d1.1 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Mahan_homozygous_nonref_SNPs_ALLSamples_pi.txt")

smpl1.1 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples.tsv")

d1.1 <- merge(d1.1, smpl1.1[, .(run, pop= type, phenotype)])
d1.1 <- merge(d1.1, geno, by = 'run')

# calculate heterozygosity at SNPs
het1.1 <- d1.1[, .(p = mean(avg_pi, na.rm=T)), by = run]

# calculate standard error
N1.1 <- d1.1[no_sites == 1, .N, by = run]

het1.1 <- merge(N1.1, het1.1)
het1.1 <- merge(het1.1, geno)

# exclude Mahan
setkey(het1.1, genotype)
lv1.1 <- rev(unlist(het1.1[, run]))
het1.1[, run := factor(run, levels = lv1.1)]



# ----- read Carya spp -----

d2 <- fread("~/workspace/heterodichogamy/Carya_spp/Carya_spp_Mahan_homozygous_nonref_SNPs_pi.txt")

# calculate heterozygosity at SNPs
het2 <- d2[, .(p = mean(avg_pi, na.rm=T)), by = pop]

# calculate standard error
N2 <- d2[no_sites == 1, .N, by = pop]

het2 <- merge(N2, het2)



# --------- Plot pecan ------


plot1 <- ggplot(het[run != "SRR15911533"], aes(x = p, y = run, color = genotype)) +
  scale_color_manual(values = c('maroon', 'tan')) +
  
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                    width = 0, 
                    position=position_dodge(width=1)) +
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = 'Heterozygosity at \npecan H-loc SNPs', y = '', color = 'Genotype', title = 'C. illinoinensis') +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(vjust = 6),
        plot.title = element_text(face = 'italic', size = 14),
        #legend.position = c(0.2,0.9)
        ) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') 
  
plot1


plot1.1 <- ggplot(het1.1[genotype != "HH" & run != 'SRR15911540'], aes(x = p, y = run, color = genotype)) +
  #scale_color_manual(values = c('maroon', 'tan')) +
  
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1)) +
  
  geom_point(data = het1.1[run == 'SRR15911540'],
             position=position_dodge(width=0.2), size = 2, color = 'black')  + 
  geom_errorbar(data = het1.1[run == 'SRR15911540'],
                aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1), color = 'black') +
  
  
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = 'Heterozygosity at \npecan H-loc SNPs', y = '', color = 'Genotype', title = 'C. illinoinensis') +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(vjust = 6),
        plot.title = element_text(face = 'italic', size = 14),
        #legend.position = c(0.2,0.9)
  ) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') 

plot1.1




# ------ Plot Carya spp -----

ENA <- c("aquatica", "cordiformis", "floridana", "glara", "laciniosa", 
         "myristiciformis", "ovata", "palmeri", "texana", "tomentosa")
EA <- c("cathayensis", "dabieshanensis", "hunanensis", 'kweichowensis', 'poilanei', "tonkinensis_1", "tonkinensis_2")

het2[pop %in% ENA, clade := 'North American']
het2[pop %in% EA, clade := 'East Asian']
het2[pop == 'glara', pop := 'glabra']
setkey(het2, clade)
lv2 <- rev(unlist(het2[, pop]))
het2[,  pop := factor(pop, levels = lv2)]

plot2 <- ggplot(het2[pop != 'suspected_ovata'], aes(x = p, y = pop, color = clade)) +
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1)) +
  scale_color_manual(values = c("gray", 'black')) +
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = 'Heterozygosity at \npecan H-loc SNPs', y = '', 
       color = 'Clade', title = 'Carya spp.') +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(b = 20, t = 5, unit = 'pt'),
        plot.title = element_text(face = 'italic', size = 14),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 5),
       # legend.position = c(0.8,0.9)
       legend.position = 'none'
        ) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') 

plot2

library(cowplot)
plot_grid(plot1, plot2, ncol = 1, align = 'v')



# ------ look at positions of SNPs -----
het_ids_NA <- het2[p > 0.2 & clade == 'North American', pop]
het_ids_EA <- het2[p > 0.2 & clade == 'East Asian', pop]


ggplot(d[pop != 'Mahan' & avg_pi == 1], aes(x = window_pos_1), y = avg_pi) + 
  geom_rug(color = 'gray', length = unit(1, 'cm')) + 
  geom_rug(data = d2[!is.na(avg_pi) & pop %in% het_ids_NA], color = 'yellow', length = unit(1.8, 'cm'))  + 
  
  geom_rug(data = d2[avg_pi == 1 & pop %in% het_ids_NA], color = 'blue', length = unit(1.1, 'cm')) + 
  geom_rug(data = d2[avg_pi == 1 & pop %in% het_ids_EA], color = 'red', length = unit(1.8, 'cm'))  + 

  geom_rug(data = d2[!is.na(avg_pi) & pop %in% het_ids_EA], color = 'green', length = unit(1.5, 'cm'))  + 
  
  theme_classic() + 
  theme(aspect.ratio = 0.2) 

d2[window_pos_1 < 6475000 & pop %in% het_ids_NA & avg_pi == 1]
d2[window_pos_1 > 6625000 & pop %in% het_ids_EA]

6474328






                