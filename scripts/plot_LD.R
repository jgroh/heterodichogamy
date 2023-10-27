library(viridis)
ld <- fread("~/workspace/heterodichogamy/regia/Jregia_founders_geno_r2.geno.ld")

#NDR1 coords 31868751	31870103
#TPP coords 31884268	31887072

setnames(ld, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))
ggplot(ld[(pos1 > 31.86e6 & pos1 < 31.895e6) &
            (pos2 > 31.86e6 & pos2 < 31.895e6)], 
       aes(x = pos1, y = pos2, fill = r2)) + 
  geom_tile(width = 150, height = 150) +
  annotate("rect", xmin = 31868751, xmax = 31870103, 
           ymin = 31868751, ymax = 31870103,
           alpha = .9,fill = 'NA', color = 'black') + 
  annotate("rect", xmin = 31884268, xmax = 31887072, 
           ymin = 31884268, ymax = 31887072,
           alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "Position (Mb)", y = "Position (Mb)", fill = expression(r^2)) +
  scale_x_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 4),
                     labels = seq(31.86, 31.89, length.out = 4)) +
  scale_y_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 4),
                     labels = seq(31.86, 31.89, length.out = 4)) +
  theme_classic() + 
  theme(aspect.ratio = 1, 
        legend.position = c(0.9, 0.2), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))
