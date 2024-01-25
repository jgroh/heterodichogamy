library(data.table)
library(ggplot2)
library(cowplot)

pecan_meta <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")
#popfile <- fread("~/workspace/heterodichogamy/pecan/WGS_data/pixy_popfile_for_heterozygosity.txt")
pecan <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/ind_het_Lakota_v1_Hloc_SNPs_pi.txt")
spp <- fread("~/workspace/heterodichogamy/Carya_spp/results/pixy/Carya_spp_pecan_SNPs_pi.txt")


pecan <- merge(pecan, pecan_meta[, .(ID, pop = variety, phenotype, genotype)])

# calculate heterozygosity at SNPs
pecan_het <- pecan[, .(p = mean(avg_pi, na.rm=T)), by = .(pop, genotype)]
spp_het <- spp[, .(p = mean(avg_pi, na.rm=T)), by = pop]


# calculate standard errors
pecan_N <- pecan[no_sites == 1, .N, by = pop]
spp_N <- spp[no_sites == 1, .N, by = pop]

pecan_het <- merge(pecan_het, pecan_N)
spp_het <- merge(spp_het, spp_N)

# ===== Plots =====

pecan_het[genotype == 'hh', genotype := 'gg']
pecan_het[genotype == 'Hh', genotype := 'Gg']
pecan_het[genotype == 'HH', genotype := 'GG']
pecan_het[, genotoype := factor(genotype, levels =c("gg", "Gh", 'GG'))]


plot1 <- ggplot(pecan_het, aes(x = p, y = pop, color = genotype)) +
  scale_color_manual(values = c('tan', 'maroon', 'turquoise4')) +
  
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1)) +
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = 'Heterozygosity at \npecan G-loc SNPs', y = '', color = '', title = 'C. illinoinensis') +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(vjust = 6),
        plot.title = element_text(face = 'italic', size = 14),
        legend.text = element_text(face = 'italic')
  ) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
  geom_vline(xintercept = 1, linetype = 2, color = 'gray') 

plot1



# ----- Plot Carya spp -----
ENA <- c("aquatica", "cordiformis", "floridana", "glara", "laciniosa", 
         "myristiciformis",  "palmeri", "texana", "tomentosa", 'ovata_1', 'ovata_2')
EA <- c("cathayensis", "dabieshanensis", "hunanensis", 'kweichowensis', 'poilanei', "tonkinensis_1", "tonkinensis_2")

spp_het[pop %in% ENA, clade := 'North American']
spp_het[pop %in% EA, clade := 'East Asian']
spp_het[pop == 'glara', pop := 'glabra']
spp_het[pop == 'tonkinensis_2', pop := 'tonkinensis']
setkey(spp_het, clade)
lv2 <- rev(unlist(spp_het[, pop]))
spp_het[,  pop := factor(pop, levels = lv2)]

plot2 <- ggplot(spp_het[!pop %in% c('unknown', 'tonkinensis_1')], aes(x = p, y = pop, color = clade)) +
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1)) +
  scale_color_manual(values = c("gray", 'black')) +
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = 'Heterozygosity at \npecan G-loc SNPs', y = '', 
       color = 'Clade', title = 'Carya spp.') +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
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
plot_grid(plot1, plot2, ncol = 2, align = 'h')


# ------ look at positions of SNPs -----
het_ids_NA <- spp_het[p > 0.2 & clade == 'North American', pop]
het_ids_EA <- spp_het[p > 0.2 & clade == 'East Asian', pop]

ggplot(spp[avg_pi == 1], aes(x = window_pos_1), y = avg_pi) + 
  #geom_rug(color = 'gray', length = unit(1, 'cm')) + 
  #geom_rug(data = spp[!is.na(avg_pi) & pop %in% het_ids_NA], color = 'yellow', length = unit(1.8, 'cm'))  + 
  #geom_rug(data = spp[!is.na(avg_pi) & pop %in% het_ids_EA], color = 'green', length = unit(1.8, 'cm'))  + 
  geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_NA], color = 'blue', length = unit(1, 'cm')) + 
  geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_EA], color = 'red', length = unit(1.1, 'cm'))  + 
  theme_classic() + 
  theme(aspect.ratio = 0.2) 

spp[window_pos_1 < 6600000 & pop %in% het_ids_NA & !is.na(avg_pi) & window_pos_1 == 6532749]

spp[window_pos_1 < 6600000 & pop %in% het_ids_EA & !is.na(avg_pi),]

6474328









