library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
require(grid)   # for the textGrob() function

# alignment coordinates
coords1 <- fread("~/workspace/heterodichogamy/01_G_locus_structure/dotplots/HJregBNU_vs_hJregCha_discontiguous_megablast.csv")        
coords2 <- fread("~/workspace/heterodichogamy/01_G_locus_structure/dotplots/HJman_vs_hJman_discontiguous_megablast.csv")
coords3 <- fread("~/workspace/heterodichogamy/01_G_locus_structure/dotplots/HJcaliPrimary_vs_hJcaliAlt_discontiguous_megablast.csv")

coords1[, Query_start := Query_start + 30760000]
coords1[, Query_end := Query_end + 30760000]

coords1[, Subject_start := Subject_start + 31868000]
coords1[, Subject_end := Subject_end + 31868000]

coords2[, Query_start := Query_start + 32670000]
coords2[, Query_end := Query_end + 32670000]

coords2[, Subject_start := Subject_start + 31760000]
coords2[, Subject_end := Subject_end + 31760000]

coords3[, Query_start := Query_start + 31358500]
coords3[, Query_end := Query_end + 31358500]

coords3[, Subject_start := Subject_start + 30498500]
coords3[, Subject_end := Subject_end + 30498500]

# locations of genes
genes_Jreg_H <- data.table(start = c(30759708, 30784000), end = c(30761060, 30786804))
genes_Jreg_h <- data.table(start = c(31868751,31884268), end = c(31870103,31887072))
  
genes_Jman_H <- data.table(start = c(32669662,32692488), end = c(32671014,32695292))
genes_Jman_h <- data.table(start = c(31760070,31776222), end = c(31761382,31779026))

genes_Jcali_H <- data.table(start = c(31358581,31380787), end = c(31359933, 31383591))
genes_Jcali_h <- data.table(start = c(30498727,30517523), end = c(30500052, 30520347))
  

# Chandler TEs
cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))

hJregTE <- fread("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/JregiaV2.fa.mod.EDTA.TEanno.gff3", skip = "NC_", col.names = cnames)
HJregTE <- fread("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/BNU/JRE_v3.3.fasta.mod.EDTA.TEanno.gff3", skip = "QKZ", col.names = cnames)


# make plots
plot1 <- ggplot(coords1) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 
  annotate("rect", xmin = genes_Jreg_H$start[1],
           xmax = genes_Jreg_H$end[1], 
           ymin=genes_Jreg_h$start[1], 
           ymax = genes_Jreg_h$end[1], 
           alpha = .3,fill = '#94435b') + 
  annotate("rect", xmin = genes_Jreg_H$start[2],
           xmax = genes_Jreg_H$end[2], 
           ymin=genes_Jreg_h$start[2], 
           ymax = genes_Jreg_h$end[2], 
           alpha = .3,fill = 'blue') +
  # 
  # annotate("rect", xmin = 30763941,
  #          xmax = 30764115, 
  #          ymin=-Inf, 
  #          ymax = Inf, 
  #          alpha = .3,fill = 'green') +
  # annotate("rect", xmin =   30765316,
  #          xmax = 30765662, 
  #          ymin=-Inf, 
  #          ymax = Inf, 
  #          alpha = .3,fill = 'green') +
  # annotate("rect", xmin =   30769251,
  #          xmax =  30769394, 
  #          ymin=-Inf, 
  #          ymax = Inf, 
  #          alpha = .3,fill = 'green') +
  # annotate("rect", xmin =   30769796,
  #          xmax = 30769965, 
  #          ymin=-Inf, 
  #          ymax = Inf, 
  #          alpha = .3,fill = 'green') +
  # annotate("rect", xmin =   30770688,
  #          xmax = 30771245, 
  #          ymin=-Inf, 
  #          ymax = Inf, 
  #          alpha = .3,fill = 'green') +
  

  
  
  theme_classic() + 
  #labs(x = 'H haplotype (Mb)', y = "h haplotype (Mb)", title = 'J. regia') +
  #labs(x = 'H haplotype (Mb)', y = "h haplotype (Mb)") +
  labs(x = "", y = "", title = 'J. regia') +
  
  scale_x_continuous(breaks = seq(30.765e6, 30.785e6, length.out = 3), labels = seq(30.765e6, 30.785e6, length.out = 3)/1e6 ) + 
  scale_y_continuous(breaks = seq(31.87e6, 31.885e6, length.out = 4), labels = sprintf("%.3f", seq(31.87e6, 31.885e6, length.out = 4)/1e6 )) +
  theme(aspect.ratio = 0.6,
      #plot.margin = unit(c(1,1,1,1), 'lines'),
      axis.text = element_text(size=10),
      #axis.title.y = element_text(size=14, vjust = 2),
      #axis.title.x = element_text(size = 14, vjust = -2),
      #text = element_text(family = "", face = 'bold'), 
      plot.title = element_text(face = 'italic', size = 10),
      #plot.margin = unit(c(.1, 0, 0, 0), "cm")
  ) + 
  #geom_rect(data = hJregTE[Chr == 'NC_049911.1' & Start > min(coords1[, Subject_start]) & End < max(coords1[, Subject_end]) ],
  #          aes(xmin = -Inf, xmax = Inf, ymin = Start, ymax = End), fill = 'lightgray') + 
  
  #geom_rect(data = HJregTE[Chr == 'chr7' & Start > min(coords1[, Query_start]) & End < max(coords1[, Query_end]) ],
  #          aes(ymin = -Inf, ymax = Inf, xmin = Start, xmax = End), fill = 'lightgray') + 
  
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 

  #geom_segment(aes(y = 31883277, yend = 31884479, x = 30785000, xend = 30785000)) + 

  # region of high Dxy used for dating estimate
  annotate("rect", xmin = -Inf,
         xmax = Inf, 
         ymin= 31871250 - 250, 
         ymax = 31871250 + 250, 
         alpha = .3,fill = 'purple') 



plot1

# make plots
plot2 <- ggplot(coords2[]) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end), 
               linewidth = 0.5) + 
  annotate("rect", xmin = genes_Jman_H$start[1],
           xmax = genes_Jman_H$end[1], 
           ymin=genes_Jman_h$start[1], 
           ymax = genes_Jman_h$end[1], 
           alpha = .3,fill = '#94435b') + 
  annotate("rect", xmin = genes_Jman_H$start[2],
           xmax = genes_Jman_H$end[2], 
           ymin=genes_Jman_h$start[2], 
           ymax = genes_Jman_h$end[2], 
           alpha = .3,fill = 'blue') +
  theme_classic() + 
  #labs(x = 'H haplotype (Mb)', y = "h haplotype (Mb)", title = 'J. regia') +
  #labs(x = 'H haplotype (Mb)', y = "h haplotype (Mb)", title = 'J. mandshurica') +
  labs(x = "", y = "", title = 'J. mandshurica') +
  
  scale_x_continuous(breaks = seq(32.67e6, 32.69e6, length.out = 3), labels = seq(32.67e6, 32.69e6, length.out = 3)/1e6 ) + 

  scale_y_continuous(breaks = seq(31.76e6, 31.775e6, length.out = 4), labels = sprintf("%.3f", seq(31.76e6, 31.775e6, length.out = 4)/1e6 )) +
  theme(aspect.ratio = 0.6,
        #plot.margin = unit(c(1,1,1,1), 'lines'),
        #plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
        axis.text = element_text(size=10),

        #axis.title.y = element_text(size=14, vjust = 2),
        #axis.title.x = element_text(size = 14, vjust = -2),
       # text = element_text(family = "", face = 'bold'), 
        plot.title = element_text(face = 'italic', size = 10)
  )

plot2




plot3 <- ggplot(coords3) + 
  geom_segment(aes(x = Query_start, xend = Query_end, 
                   y = Subject_start, yend = Subject_end),
               linewidth = 0.5) + 
  annotate("rect", xmin = genes_Jcali_H$start[1],
           xmax = genes_Jcali_H$end[1], 
           ymin=genes_Jcali_h$start[1], 
           ymax = genes_Jcali_h$end[1], 
           alpha = .3,fill = '#94435b') + 
  annotate("rect", xmin = genes_Jcali_H$start[2],
           xmax = genes_Jcali_H$end[2], 
           ymin=genes_Jcali_h$start[2], 
           ymax = genes_Jcali_h$end[2], 
           alpha = .3,fill = 'blue') + 
  annotate("rect", xmin = 31361145,
           xmax = 31361188, 
           ymin=genes_Jcali_h$start[2], 
           ymax = genes_Jcali_h$end[2], 
           alpha = .3,fill = 'orange') +
  annotate("rect", xmin = 31369297,
           xmax = 31369348, 
           ymin=genes_Jcali_h$start[2], 
           ymax = genes_Jcali_h$end[2], 
           alpha = .3,fill = 'orange') +
  #geom_vline(aes(xintercept = 31369300)) + 
  #geom_vline(aes(xintercept = 31381130)) + 
  #geom_vline(aes(xintercept = 31382320)) + 
  #geom_vline(aes(xintercept = 31383130)) + 
  
  
  theme_classic() + 
  #labs(x = 'H Jcali (Mb)', y = "h Jcali (Mb)") +
  #labs(x = 'H haplotype (Mb)', y = "h haplotype (Mb)", title = 'J. californica') +
  labs(x = "", y = "", title = 'J. californica') +
  
  scale_x_continuous(breaks = seq(31.36e6, 31.38e6, length.out = 3), labels = seq(31.36e6, 31.38e6, length.out = 3)/1e6 ) + 
  scale_y_continuous(breaks = seq(30.5e6, 30.52e6, length.out = 3), labels = sprintf("%.2f",seq(30.5e6, 30.52e6, length.out = 3)/1e6 )) +
  theme(aspect.ratio = 0.6,
        #plot.margin = unit(c(1,1,1,1), 'lines'),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        axis.text = element_text(size=10),
        #axis.title.y = element_text(size=14, vjust = 2),
        #axis.title.x = element_text(size = 14, vjust = -2),
        #text = element_text(family = "", face = 'bold'),
        plot.title = element_text(face = "italic",size = 10)
  )

plot3

ggarrange(plot1, plot2, plot3, nrow=3, ncol=1)

figure_v <- ggarrange(plot1, plot2, plot3,
          ncol = 1, nrow = 3)
figure_v

figure_h <- ggarrange(plot1, plot2, plot3,
                      ncol = 3, nrow = 1)
figure_h

ggarrange(plot1, plot2, plot3)


  

  
