library(viridis)
library(data.table)
library(ggplot2)
library(cowplot)
regia <- fread("~/workspace/heterodichogamy/regia/Jregia_founders_geno_r2.geno.ld")
hindsii <- fread("~/workspace/heterodichogamy/hindsii/Putah_Hloc.geno.ld")

setnames(hindsii, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))

#31380614, 31383466
#31358581, 31360135
hindsii_plt <- ggplot(hindsii[(pos1 > 31.35e6 & pos1 < 31.39e6) &
               (pos2 > 31.35e6 & pos2 < 31.39e6)], 
       aes(x = pos2, y = pos1, fill = r2)) + 
  geom_tile(width = 150, height = 150) +
  annotate("rect", xmin = 31380614, xmax = 31383466, 
           ymin = 31380614, ymax = 31383466,
           alpha = .9,fill = 'NA', color = 'black') + 
  annotate("rect", xmin = 31358581, xmax = 31360135, 
           ymin = 31358581, ymax = 31360135,
           alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "Position (Mb)", y = "Position (Mb)", fill = expression(r^2), title = 'J. hindsii') +
  scale_x_continuous(breaks = seq(31.35e6, 31.39e6, length.out = 3),
                     labels = seq(31.35, 31.39, length.out = 3)) +
  scale_y_continuous(breaks = seq(31.35e6, 31.39e6, length.out = 3),
                     labels = seq(31.35, 31.39, length.out = 3)) +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        legend.position = c(0.1, 0.7), 
       # legend.key.size = unit(0.3, 'cm'),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
       plot.title = element_text(face = 'italic'))

hindsii

# ---- hindsii, just across TPP
ggplot(hindsii[(pos1 > 31.38e6 & pos1 < 31.39e6) &
                 (pos2 > 31.38e6 & pos2 < 31.39e6)], 
       aes(x = pos2, y = pos1, fill = r2)) + 
  geom_tile(width = 150, height = 150) +
  annotate("rect", xmin = 31380614, xmax = 31383466, 
           ymin = 31380614, ymax = 31383466,
           alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "Position (Mb)", y = "Position (Mb)", fill = expression(r^2)) +
  scale_x_continuous(breaks = seq(31.35e6, 31.39e6, length.out = 3),
                     labels = seq(31.35, 31.39, length.out = 3)) +
  scale_y_continuous(breaks = seq(31.35e6, 31.39e6, length.out = 3),
                     labels = seq(31.35, 31.39, length.out = 3)) +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        legend.position = c(0.1, 0.7), 
        legend.key.size = unit(0.3, 'cm'),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12))










# ------regia, plot LD mapped to Chandler ------
#NDR1 coords 31868751	31870103
#TPP coords 31884268	31887072
setnames(regia, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))
regia_plt <- ggplot(regia[(pos1 > 31.86e6 & pos1 < 31.895e6) &
            (pos2 > 31.86e6 & pos2 < 31.895e6)], 
       aes(x = pos2, y = pos1, fill = r2)) + 
  geom_tile(width = 150, height = 150) +
  annotate("rect", xmin = 31868751, xmax = 31870103, 
           ymin = 31868751, ymax = 31870103,
           alpha = .9,fill = 'NA', color = 'black') + 
  annotate("rect", xmin = 31884268, xmax = 31887072, 
           ymin = 31884268, ymax = 31887072,
           alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "Position (Mb)", y = "Position (Mb)", fill = expression(r^2), title = 'J. regia') +
  scale_x_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 4),
                     labels = seq(31.86, 31.89, length.out = 4)) +
  scale_y_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 4),
                     labels = seq(31.86, 31.89, length.out = 4)) +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        legend.position = c(0.1, 0.7), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(face = 'italic'))


plot_grid(hindsii_plt, regia_plt)
