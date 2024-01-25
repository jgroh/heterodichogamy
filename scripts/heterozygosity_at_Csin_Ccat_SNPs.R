library(data.table)
library(ggplot2)
library(cowplot)

# ---- genes -----
genes <- fread("~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_Csin/Chr4_gene_coords.txt")
# define variable to indicate whether gene is in locus or outside
genes[, G_loc := ifelse(V4 >= 6074902 & V5 <= 6421582, 'G', 'BG')]
setkey(genes, V4)

#pecan_meta <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")
#popfile <- fread("~/workspace/heterodichogamy/pecan/WGS_data/pixy_popfile_for_heterozygosity.txt")
#pecan <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/ind_het_Lakota_v1_Hloc_SNPs_pi.txt")

#spp <- fread("~/workspace/heterodichogamy/Carya_spp/Carya_spp_Csin_Ccat_chr4_SNPs_pi.txt")
spp <- fread("~/workspace/heterodichogamy/Carya_spp/Carya_spp_Csin_Ccat_chr4_SNPs_DP6_GQ20_pi.txt")

spp_Gloc <- spp[window_pos_1 > 6074902 & window_pos_2 < 6421582]
spp_3LHS <- spp[window_pos_1 > 6074900 & window_pos_2 < 6088683]
spp_1RHS <- spp[window_pos_1 > 6419954 & window_pos_2 < 6421582]


#pecan <- merge(pecan, pecan_meta[, .(ID, pop = variety, phenotype, genotype)])

# calculate heterozygosity at SNPs
#pecan_het <- pecan[, .(p = mean(avg_pi, na.rm=T)), by = .(pop, genotype)]
spp_het_Gloc <- spp_Gloc[, .(p = mean(avg_pi, na.rm=T)), by = pop]
spp_het_chr4 <- spp[, .(p = mean(avg_pi, na.rm=T)), by = pop]
spp_het_3LHS <- spp_3LHS[, .(p = mean(avg_pi, na.rm=T)), by = pop]
spp_het_1RHS <- spp_1RHS[, .(p = mean(avg_pi, na.rm=T)), by = pop]




# calculate standard errors
#pecan_N <- pecan[no_sites == 1, .N, by = pop]

spp_N_chr4 <- spp[no_sites == 1, .N, by = pop]
spp_N_Gloc <- spp_Gloc[no_sites == 1, .N, by = pop]
spp_N_3LHS <- spp_3LHS[no_sites == 1, .N, by = pop]
spp_N_1RHS <- spp_1RHS[no_sites == 1, .N, by = pop]


#pecan_het <- merge(pecan_het, pecan_N)

spp_het_chr4 <- merge(spp_het_chr4, spp_N_chr4)
spp_het_Gloc <- merge(spp_het_Gloc, spp_N_Gloc)
spp_het_3LHS <- merge(spp_het_3LHS, spp_N_3LHS)
spp_het_1RHS <- merge(spp_het_1RHS, spp_N_1RHS)

# ===== Plots =====

# pecan_het[genotype == 'hh', genotype := 'gg']
# pecan_het[genotype == 'Hh', genotype := 'Gg']
# pecan_het[genotype == 'HH', genotype := 'GG']
# pecan_het[, genotoype := factor(genotype, levels =c("gg", "Gh", 'GG'))]


# plot1 <- ggplot(pecan_het, aes(x = p, y = pop, color = genotype)) +
#   scale_color_manual(values = c('tan', 'maroon', 'turquoise4')) +
#   
#   geom_point(position=position_dodge(width=0.2), size = 2)  + 
#   geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
#                     xmax = p + 1.96*sqrt(p*(1-p)/N)), 
#                 width = 0, 
#                 position=position_dodge(width=1)) +
#   scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
#   labs(x = 'Heterozygosity at \npecan G-loc SNPs', y = '', color = '', title = 'C. illinoinensis') +
#   theme_classic() + 
#   theme(aspect.ratio = 1, 
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title.y = element_text(vjust = 6),
#         plot.title = element_text(face = 'italic', size = 14),
#         legend.text = element_text(face = 'italic')
#   ) + 
#   geom_vline(xintercept = 0, linetype = 2, color = 'gray') + 
#   geom_vline(xintercept = 1, linetype = 2, color = 'gray') 
# 
# plot1



# ----- Plot Carya spp -----
ENA <- c("aquatica", "cordiformis", "floridana", "glara", "laciniosa", 
         "myristiciformis",  "palmeri", "texana", "tomentosa", 'ovata_1', 'ovata_2')
EA <- c("cathayensis", "dabieshanensis", "hunanensis", 'kweichowensis', 'poilanei', "tonkinensis_1", "tonkinensis_2")

spp_het_chr4[pop %in% ENA, clade := 'North American']
spp_het_chr4[pop %in% EA, clade := 'East Asian']
spp_het_chr4[pop == 'glara', pop := 'glabra']
spp_het_chr4[pop == 'tonkinensis_2', pop := 'tonkinensis']

setkey(spp_het_chr4, clade)

spp_het_Gloc[pop %in% ENA, clade := 'North American']
spp_het_Gloc[pop %in% EA, clade := 'East Asian']
spp_het_Gloc[pop == 'glara', pop := 'glabra']
spp_het_Gloc[pop == 'tonkinensis_2', pop := 'tonkinensis']

setkey(spp_het_Gloc, clade)

spp_het_3LHS[pop %in% ENA, clade := 'North American']
spp_het_3LHS[pop %in% EA, clade := 'East Asian']
spp_het_3LHS[pop == 'glara', pop := 'glabra']
spp_het_3LHS[pop == 'tonkinensis_2', pop := 'tonkinensis']

setkey(spp_het_3LHS, clade)

spp_het_1RHS[pop %in% ENA, clade := 'North American']
spp_het_1RHS[pop %in% EA, clade := 'East Asian']
spp_het_1RHS[pop == 'glara', pop := 'glabra']
spp_het_1RHS[pop == 'tonkinensis_2', pop := 'tonkinensis']

setkey(spp_het_1RHS, clade)

lv2 <- rev(unlist(spp_het_Gloc[, pop]))
spp_het_Gloc[,  pop := factor(pop, levels = lv2)]
spp_het_3LHS[,  pop := factor(pop, levels = lv2)]
spp_het_1RHS[,  pop := factor(pop, levels = lv2)]


# all chr 4
ggplot(spp_het_chr4[!pop %in% c('unknown', 'tonkinensis_1')], aes(x = p, y = pop, color = clade)) +
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

# G-locus

ggplot(spp_het_Gloc[!pop %in% c('unknown', 'tonkinensis_1', 'dabieshanensis')], aes(x = p, y = pop, color = clade)) +
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

# 3 genes at LHS
ggplot(spp_het_3LHS[!pop %in% c('unknown', 'tonkinensis_1', 'dabieshanensis')], aes(x = p, y = pop, color = clade)) +
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


# gene at RHS
ggplot(spp_het_1RHS[!pop %in% c('unknown', 'tonkinensis_1', 'dabieshanensis')], aes(x = p, y = pop, color = clade)) +
  geom_point(position=position_dodge(width=0.2), size = 2)  + 
  geom_errorbar(aes(xmin = p - 1.96*sqrt(p*(1-p)/N),
                    xmax = p + 1.96*sqrt(p*(1-p)/N)), 
                width = 0, 
                position=position_dodge(width=1)) +
  scale_color_manual(values = c("gray", 'black')) +
  scale_x_continuous(breaks = c(0, .25, 0.5,0.75, 1)) +
  labs(x = '', y = '', 
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






# plot along chromosome
# Assign each bp to its containing window 
round_up <- function(x) {
  result <- ceiling(x / window_size) * window_size
  return(result)
}

round_down <- function(x) {
  result <- floor(x / window_size) * window_size
  return(result)
}

window_size <- 5e3
spp[, window := cut(window_pos_1, 
                      breaks = seq(from = round_down(window_pos_1[1]), to = round_up(max(window_pos_1)), by = window_size) , 
                      include.lowest=T, 
                      labels = seq(from = round_down(window_pos_1[1]), 
                                   to = round_up(max(window_pos_1)) - window_size, by = window_size) + 0.5*window_size), 
      by = pop]


window_pi <- spp[, .(window_pi = mean(avg_pi, na.rm=T)), by = .(pop,window)]
window_pi[, window := as.numeric(as.character(window))]

ggplot(window_pi[pop == 'kweichowensis' & window > 5.6e6 & window < 6.9e6], aes(x = window, y = window_pi)) + geom_point()




# ========= kweichowensis =====
# get cds for each gene 
k <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'kweichowensis']
aquatica <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'aquatica']
poi <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'poilanei']
hun <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'hunanensis']
cat <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'cathayensis']
ov1 <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'ovata_1']
ov2 <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'ovata_2']
tk <- spp[chromosome == 'Chr04' & window_pos_1 > 5000000 & window_pos_2 < 7000000 & pop == 'tonkinensis']

focal_genes <- genes[V4 > 5000000 & V5 < 7000000]

for(i in 1:nrow(focal_genes)){
  start <- focal_genes[i, V4]
  end <- focal_genes[i, V5]
  
  k_val <- k[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  a_val <- aquatica[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  poi_val <- poi[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  hun_val <- hun[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  cat_val <- cat[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  ov1_val <- ov1[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  ov2_val <- ov2[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
  tk_val <- tk[window_pos_1 >= start & window_pos_1 <= end, mean(avg_pi, na.rm=T)]
    
    
  focal_genes[i, k_pi := k_val]
  focal_genes[i, aquatica_pi := a_val]
  focal_genes[i, poi_pi := poi_val]
  focal_genes[i, hun_pi := hun_val]
  focal_genes[i, cat_pi := cat_val]
  focal_genes[i, ov1_pi := ov1_val]
  focal_genes[i, ov2_pi := ov2_val]
  focal_genes[i, tk_pi := tk_val]
  
}
#focal_genes[, k_pi]



kweichowensis_plt <- ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = k_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = V4, y = k_pi)) + 
  labs(y = '', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. kweichowensis') + # if bottom
  scale_x_continuous(breaks = seq(5.7e6, 6.7e6, length.out = 6), labels = seq(5.7, 6.7, length.out = 6) )  +
  
  theme_classic() + 
  theme(aspect.ratio = .2,
        #plot.margin = margin(0,0,0,0, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(face = 'italic'),
        #legend.position = c(0.15, 0.85),
        legend.position = 'none',
        legend.key.size = unit(0.3, 'cm'),
  )   




aquatica_plt <- ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = aquatica_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = (V4 + V5)/2, y = aquatica_pi)) + 
  labs(y = '', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. aquatica') + # if bottom
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


# ovata 1
ovata1 <- ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = ov1_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  labs(y = '', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. ovata') + # if bottom
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

ovata1


# ovata 2
ovata2 <- ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = ov2_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  labs(y = '', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. ovata') + # if bottom
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

ovata2

# cathayensis
cathayensis_plt <- ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = cat_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = (V4 + V5)/2, y = cat_pi)) + 
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


plot_grid(kweichowensis_plt, cathayensis_plt, ovata1, ovata2, ncol = 2)



# ========= 3 genes LHS  =======
# Coords of Protein phosphate inhibitor 2 like (ASI021694) Chr04:6074901-6077206
# Coords of LOC122306891_uncharacterized (ASI021582) Chr04:6080732-6082069
# Coords of EMS1-like (ASI021582) Chr04:6084797-6088682

# look at heterozygosity for the three LHS SNPs in th
spp_3LHS


spp[chromosome == 'Chr04' & window_pos_1 > 6074900 & window_pos_2 < 6088683]

focal_genes[V4 > 6074900 & V5 < 6088683]





ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = k_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = V4, y = k_pi)) + 
  labs(y = 'Heterozygosity\nat Csin-Ccat SNPs', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. kweichowensis') + # if bottom
  scale_x_continuous(breaks = seq(5.7e6, 6.7e6, length.out = 6), labels = seq(5.7, 6.7, length.out = 6) )  +
  
  theme_classic() + 
  theme(aspect.ratio = .2,
        #plot.margin = margin(0,0,0,0, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(face = 'italic'),
        #legend.position = c(0.15, 0.85),
        legend.position = 'none',
        legend.key.size = unit(0.3, 'cm'),
  )   




# ------------ Others --------


# cathayensis
ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = cat_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = (V4 + V5)/2, y = cat_pi)) + 
  labs(y = 'Heterozygosity\nat Csin-Ccat SNPs', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. cathayensis') + # if bottom
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


# poilanei
ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = poi_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = (V4 + V5)/2, y = poi_pi)) + 
  labs(y = 'Heterozygosity\nat Csin-Ccat SNPs', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)'), title = 'C. poilanei') + # if bottom
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


# hunanensis
ggplot(focal_genes[V4 > 5700000 & V4 < 6900000]) + 
  geom_point(aes(x = (V4 + V5)/2, y = hun_pi)) + 
  scale_y_continuous(limits = c(-.1, 1)) + 
  geom_rect(data = genes[V4 > 5700000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0, fill = G_loc)) + 
  #geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) + 
  geom_point(aes(x = (V4 + V5)/2, y = hun_pi)) + 
  labs(y = 'Heterozygosity\nat Csin-Ccat SNPs', x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)')) + # if bottom
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





