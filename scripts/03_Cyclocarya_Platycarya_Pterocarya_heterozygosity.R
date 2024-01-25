library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(readxl)

# ----- Cylocarya ----- 
Cyc_pi_HC <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/Cpal_diploids_chr4_6-7.2Mb_pi.txt")

# ---- Carya ----
Car_pi_HC <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars_chr4_1kb_pi.txt")
smpl <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_reassigned.tsv")
geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c("run", 'genotype'))

Car_pi_HC <- merge(Car_pi_HC, smpl[, .(run, pop= type, phenotype)])
Car_pi_HC <- merge(Car_pi_HC, geno, by = 'run')
Car_pi_HC[, pos := (window_pos_1 + window_pos_2)/2]
Car_pi_HC[genotype == 'HH', genotype := 'GG']
Car_pi_HC[genotype == 'Hh', genotype := 'Gg']
Car_pi_HC[genotype == 'hh', genotype := 'gg']

# # ---- all Carya samples
# Car_pi_HC2 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars_chr4_1kb_ALL_pi.txt")
# smpl2 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples.tsv")
# 
# Car_pi_HC2 <- merge(Car_pi_HC2, smpl2[, .(run, pop= type, phenotype)])
# Car_pi_HC2 <- merge(Car_pi_HC2, geno, by = 'run')
# Car_pi_HC2[, pos := (window_pos_1 + window_pos_2)/2]
# 

# ----- Platycarya -----
Pstr <- fread("~/workspace/heterodichogamy/Platycarya_strobilaceae/Pstr_Chr04_pi.txt")

# ----- Pterocarya -----
Pste_pi_HC <- fread("~/workspace/heterodichogamy/Pterocarya_stenoptera/results/pixy/Pste_Carya_loc_1kb_pi.txt")




# ----- Carya locus ----- 




HC1 <- ggplot(Car_pi_HC[window_pos_1 > 6380000 & window_pos_2 < 6740000], 
              aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = genotype)) + 
  geom_point(alpha = 0.5, size = 0.8) + 
  scale_color_manual(values = c( 'tan', 'maroon', 'turquoise4')) +
  geom_smooth(aes(group = pop), method = 'loess', se = F, span = 0.4, linewidth = 0.8) + 
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Pecan') +
  geom_vline(xintercept = 6461404, linetype = 2) + 
  geom_vline(xintercept = 6658636, linetype = 2) + 
  
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 60, hjust = 0.8, vjust = 0.8),
        legend.text = element_text(face = 'italic'),
        plot.title=element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
HC1

HC2 <- ggplot(Cyc_pi_HC[window_pos_1 > 6430000 & window_pos_2 < 6840000], 
       aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  geom_point(alpha = 0.5, size =0.8) + 
  geom_smooth(method = 'loess', se = F, span = 0.4, linewidth = 0.8) + 
  scale_color_brewer(type = 'qual', palette = 'Paired') + 
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Cyclocarya paliurus') +
  geom_vline(xintercept = 6747937, linetype = 2) + 
  geom_vline(xintercept = 6553915, linetype = 2) + 

  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 60, hjust = 0.8, vjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 

HC2

HC3 <- ggplot(Pste_pi_HC[window_pos_1 > 10000 & window_pos_2 < 500000], 
              aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  geom_point(alpha = 0.5, size =0.8) + 
  geom_smooth(method = 'loess', se = F, span = 0.4, linewidth = 0.8) + 
  #scale_color_brewer(type = 'qual', palette = 'Paired') + 
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Pterocarya stenoptera') +
  geom_vline(xintercept = 148841, linetype = 2) + 
  geom_vline(xintercept = 322258, linetype = 2) + 
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 60, hjust = 0.8, vjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
HC3


HC4 <- ggplot(Pstr[window_pos_1 > 26000000 & window_pos_2 < 28000000], 
       aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  geom_point(alpha = 0.4, size =0.8) + 
  geom_smooth(method = 'loess', se = F, span = 0.4, linewidth = 0.8) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') + 
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Platycarya strobilaceae') +
  geom_vline(xintercept = 27326709, linetype = 2) + 
  geom_vline(xintercept = 26558687, linetype = 2) + 
  
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 60, hjust = 0.8, vjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
HC4



ggarrange(HC1, HC2, HC3, HC4, align = 'v')




# ----- compare with Juglans locus -----
# ---- cylocarya
Cyc_pi_H <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/Cpal_diploids_TPPDchrom_pi.txt")

# ---- juglans regia
Jreg_pi_H <- fread("~/workspace/heterodichogamy/regia/founders_chr11_30-32Mb_500bp_pi.txt")
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", header = F, col.names = c("pop", 'phenotype', 'species'))
Jreg_pi_H <- merge(Jreg_pi_H, pheno)
Jreg_pi_H[phenotype == 'protandrous', genotype := 'gg']
Jreg_pi_H[phenotype == 'protogynous', genotype := 'Gg']
Jreg_pi_H[pop == 'JG0026', genotype := 'GG']

# ---- juglans hindsii
hindsii_pi <- fread("~/workspace/heterodichogamy/hindsii/chr11_500bp_pi.txt.gz")
hindsii_pi <- merge(hindsii_pi, pheno)
hindsii_pi[phenotype == 'protandrous', genotype := 'gg']
hindsii_pi[phenotype == 'protogynous', genotype := 'Gg']

hindsii_100kb_pi <- fread("~/workspace/heterodichogamy/hindsii/chr11_100kb_pi.txt.gz")
hindsii_100kb_pi <- merge(hindsii_100kb_pi, pheno)
hindsii_100kb_pi[phenotype == 'protandrous', genotype := 'gg']
hindsii_100kb_pi[phenotype == 'protogynous', genotype := 'Gg']

# ----- Pterocarya
Pste_pi_H <- fread("~/workspace/heterodichogamy/Pterocarya_stenoptera/results/pixy/Pste_Juglans_loc_500bp_pi.txt")


# ---- Platycarya
Pstr_pi_H <- fread("~/workspace/heterodichogamy/Platycarya_strobilaceae/Pstr_Juglans_locus_pi.txt")

# #31840000
# H1 <- ggplot(Jreg_pi_H[window_pos_1 > 31800000 & window_pos_2 < 31920000], 
#              aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = genotype)) + 
#   annotate("rect", xmin = 31884268, xmax = 31887072, ymin = 0, ymax = 1,
#            alpha = .2,fill = 'blue') +
#   annotate("rect", xmin = 31868753, xmax = 31870101, ymin = 0, ymax = 1,
#            alpha = .2,fill = 'maroon') +
#   geom_point(size = 0.8) + 
#   geom_smooth(method = 'loess', se = F, aes(group = pop), linewidth = 0.8) + 
#   scale_color_manual(values = c("tan", "turquoise4", 'maroon')) +
#   labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'J. regia') +
# 
#   theme_classic() + 
#   theme(aspect.ratio = 1, 
#         plot.title = element_text(face = 'italic', size = 15),
#         legend.text = element_text(face = 'italic'),
#         axis.text = element_text(size = 10),
#         axis.text.x = element_text(angle = 90),
#         axis.title = element_text(size = 12)) 
# H1

# look at hindsii pi genome-wide 
#ggplot(hindsii_100kb_pi, aes(x = pop, y = avg_pi)) + geom_boxplot() + 
#  theme(axis.text.x = element_text(angle = 90))


H1 <- ggplot(hindsii_pi[window_pos_1 > 31300000 & window_pos_2 < 31450000 & !is.na(genotype)],
             aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = genotype)) +
  annotate("rect", xmin = 31380787, xmax = 31383637, ymin = 0, ymax = 1,
           alpha = .2,fill = 'blue') +
  annotate("rect", xmin = 31358581, xmax = 31360135, ymin = 0, ymax = 1,
           alpha = .2,fill = 'maroon') +
  geom_point(size = 0.8) +
  geom_smooth(method = 'loess', se = F, aes(group = pop), linewidth = 0.8, span = 0.5) +
  scale_color_manual(values = c("tan", 'maroon')) +
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Juglans hindsii') +

  theme_classic() +
  theme(aspect.ratio = 1,
        plot.title = element_text(face = 'italic', size = 12),
        legend.text = element_text(face = 'italic'),
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.8, hjust = 0.8),
        axis.title = element_text(size = 12))
H1


H2 <- ggplot(Cyc_pi_H[window_pos_1 > 5500000 & window_pos_2 < 5650000], 
       aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  annotate("rect", xmin = 5570161, xmax = 5572868, ymin = 0, ymax = .15,
           alpha = .2,fill = 'blue') +
  annotate("rect", xmin = 5581882, xmax = 5580884, ymin = 0, ymax = .15,
           alpha = .2,fill = 'maroon') +
  geom_point(size = 0.8) + 
  scale_x_reverse() +
  geom_smooth(method = 'loess', se = F, linewidth = 0.8) + 
  scale_color_brewer(type = 'qual', palette = 'Paired') + 
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Cyclocarya paliurus') +
  
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.8, hjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 

H2

H3 <- ggplot(Pste_pi_H[window_pos_1 > 580000 & window_pos_2 < 660000], 
             aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  annotate("rect", xmin = 610734, xmax = 613211, ymin = 0, ymax = .4,
           alpha = .2,fill = 'blue') +
  annotate("rect", xmin = 621258, xmax = 622681, ymin = 0, ymax = .4,
           alpha = .2,fill = 'maroon') +
  geom_point(size = 0.8) + 
  scale_x_reverse() +
  geom_smooth(method = 'loess', se = F, linewidth = 0.8) + 
  #scale_color_brewer(type = 'qual', palette = 'Paired') + 
  #scale_color_manual(pal = rainbow(13)) +
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Pterocarya stenoptera') +
  
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.8, hjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 

H3


H4 <- ggplot(Pstr_pi_H[window_pos_1 > 37850000 & window_pos_2 < 38000000], 
             aes(x = (window_pos_1 + window_pos_2)/2, y = avg_pi, color = pop)) + 
  annotate("rect", xmin = 37931099, xmax = 37933477, ymin = 0, ymax = .3,
           alpha = .2,fill = 'blue') +
  annotate("rect", xmin = 37913007, xmax = 37913825, ymin = 0, ymax = .3,
           alpha = .2,fill = 'maroon') +
  geom_point(size = 0.8) + 
  geom_smooth(method = 'loess', se = F, linewidth = 0.8) + 
  #scale_color_brewer(type = 'qual', palette = 'Paired') + 
  #scale_color_manual(pal = rainbow(13)) +
  labs(x = 'Position', y = 'Heterozygosity', color = '', title = 'Platycarya strobilaceae') +
  
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.8, hjust = 0.8),
        plot.title = element_text(face = 'italic', size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 

H4




plot_grid(H1, H2, H3, H4, align = 'v')


