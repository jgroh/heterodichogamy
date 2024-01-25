library(data.table)

# ---- genes -----
genes <- fread("~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_Csin/Chr4_gene_coords.txt")
# define variable to indicate whether gene is in locus or outside
genes[, G_loc := ifelse(V4 >= 6074902 & V5 <= 6421582, 'G', 'BG')]


# ----- divergence between Carya sinensis and Carya cathayensis assemblies ------
dxy1 <- fread("~/workspace/heterodichogamy/whole_genome_alignments/QRY_vs_Csin/Csin_Ccat_dxy.txt")


# mean and outliers
dxy1[, avgDxy := weighted.mean(avg_dxy, w = window_pos_2 - window_pos_1, na.rm=T)]
dxy1[, qnt0.95 := quantile(avg_dxy, probs = c(0.95), na.rm=T)]
dxy1[, qnt0.99 := quantile(avg_dxy, probs = c(0.99), na.rm=T)]
dxy1[, outlier := ifelse(avg_dxy > qnt0.95, 1, 0)]

dxy1[chromosome == 'Chr04']

#6074902 & End < 6421582
dxy1[chromosome == 'Chr04' & window_pos_1 > 6074902 & window_pos_2 < 6421582, weighted.mean(avg_dxy, w = window_pos_2 - window_pos_1, na.rm = T)]



ggplot(dxy1[chromosome == 'Chr04' & window_pos_1 > 5600000 & window_pos_2 < 6900000]) + 
  
  #geom_rect(aes(xmin=-Inf, xmax= Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_hline(aes(yintercept = avgDxy), linetype = 2, color = 'darkgray') +
  geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  
  geom_smooth(method = 'loess', color = 'red', span = 0.4, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy, weight = abs(window_pos_2 - window_pos_1))) + 
  #scale_x_continuous(breaks = seq(6.4e6, 7e6, length.out = 4), labels = seq(6.4, 7, length.out = 4))  + 
  scale_x_continuous(breaks = seq(5.7e6, 6.7e6, length.out = 6), labels = seq(5.7, 6.7, length.out = 6) )  +
  theme_classic() + 
  scale_y_continuous(limits = c(-.01, 0.08)) + 
  labs(y = expression('D'[xy]), x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)')) + # if bottom
  #geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_NA], aes(x = window_pos_1, y = NULL), color = 'blue', length = unit(.3, 'cm')) + 
  #geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_EA], aes(x = window_pos_1, y = NULL), color = 'red', length = unit(.3, 'cm')) + 
  #geom_rug(data = spp[avg_pi == 0 & pop %in% het_ids_EA], aes(x = window_pos_1, y = NULL), color = 'green', length = unit(.3, 'cm')) + 
  geom_vline(aes(xintercept = 6074902), linetype = 2) + 
  geom_vline(aes(xintercept = 6421582), linetype = 2) + 
  
  theme(aspect.ratio = .2,
        #plot.margin = margin(0,0,0,0, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        #legend.position = c(0.15, 0.85),
        legend.position = 'none',
        legend.key.size = unit(0.3, 'cm'),
  )   + 
  geom_rect(data = genes[V4 > 5600000 & V5 < 6900000], aes(xmin = V4, xmax = V5, ymin = -0.01, ymax = 0, fill = G_loc)) + 
  geom_point(size = 0.8, aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) +  
  scale_fill_manual(values = c("gray", 'red')) 
