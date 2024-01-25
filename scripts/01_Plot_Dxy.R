library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggpubr)
require(grid)

# ===== Read Data =====
alignments_dir <- "/Users/Jeff/workspace/heterodichogamy/whole_genome_alignments"

TEs <- fread("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/JregiaV2.fa.mod.EDTA.TEanno.gff3", 
             header = F, skip = "#")



dxy <- rbindlist(lapply(list.files(alignments_dir,
                  recursive = TRUE,
                  pattern = "^Dxy.*\\.txt$",
                  full.names = TRUE), 
       fread))


dxy[, hap := substr(qryGnom, 1, 1)]
dxy[, species := ifelse(substr(qryGnom,1,2) %in% c('hJ', 'HJ'), 
             substr(qryGnom,2,5), substr(qryGnom,1,4))]
dxy[substr(species,1,1) != 'J', hap := 'N']


# ===== Calculate Mean and Quantiles of Dxy =====
dxy[, avgDxy := mean(Dxy, na.rm=T), by = .(refGnom, qryGnom)]
dxy[, qnt0.95 := quantile(Dxy, probs = c(0.95), na.rm=T), by = .(refGnom, qryGnom)]
dxy[, qnt0.99 := quantile(Dxy, probs = c(0.99), na.rm=T), by = .(refGnom, qryGnom)]

# ===== Plot divergence of all alignments against Chandler =====

DxyhJregSub <- dxy[refGnom == 'hJregCha' & 
                   
                     window > 31.855e6 & 
                     window < 31.9e6]
DxyhJregSub[species == 'Cycl', species := 'Cpal']
DxyhJregSub[, species1 := factor(species, levels = c("Jreg", 'Jsig', 'Jman', 'Jcal', 'Jmic', 'Jnig', 'Pste', 'Cpal', 'Pstr', 'Plon', 'Cill'))]

DxyhJregSub[hap == 'h', hap := 'g']
DxyhJregSub[hap == 'H', hap := 'G']

ggplot(DxyhJregSub, aes(x = window, y = Dxy, color = hap, group = qryGnom) ) + 
  facet_wrap(~species1) +
  #geom_point()  +
  labs(x = "g haplotype (Mb)", y = 'Nucleotide divergence', color = '') +
  geom_rect(data = DxyhJregSub, aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.95), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line() +
  
  scale_x_continuous(breaks = seq(31.868e6, 31.888e6, length.out = 3), labels = seq(31.868e6, 31.888e6, length.out = 3)/1e6 ) +
  scale_color_manual(values = c('darkgoldenrod1', "#413a6e", 'darkgray')) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
        legend.position = c(0.9, 0.1), 
        legend.key.size = unit(.5, 'cm'), 
        legend.text = element_text(size = 14), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 14)) + 
  # NDR1
   annotate("rect", xmin = 31868750,
            xmax = 31869735, 
            ymin=0, 
            ymax = 0.3, 
            alpha = .3,fill = 'turquoise4')  +
  #TPPD
  annotate("rect", xmin = 31884475,
           xmax = 31887080, 
           ymin=0, 
           ymax = 0.3, 
           alpha = .3,fill = '#94435b') + 
  #HJ3
  #annotate("rect", xmin = 31.868e6+3067,
   #        xmax = 31.868e6 + 3213, 
   #        ymin=0, 
   #        ymax = 0.3,fill = 'tan') +
  # mean
  geom_hline(data = DxyhJregSub, aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) 



#

# ===== Plot locations of Annotated transposons =====

Dxy_TEs <- dxy[refGnom == 'hJregCha' & 
                # species == 'Jreg' &
                     window > 31.86e6 & 
                     window < 31.898e6]

TEs_sub <- TEs[V1 == 'NC_049911.1' & V4 > 31.868e6 & V5 < 31.898e6]
TEs_sub[, N := seq_len(.N)/.N]

ggplot(Dxy_TEs[hap == 'h'] ) + 

  #geom_point()  +
  #facet_wrap(~species) +
  labs(x = "Position", y = 'Divergence against Chandler (Jreg h)') +
  geom_line(aes(x = window, y = Dxy, color = hap, group = qryGnom)) +
  
  scale_x_continuous(breaks = seq(31.868e6, 31.888e6, length.out = 3), labels = seq(31.868e6, 31.888e6, length.out = 3)/1e6 ) +
  scale_color_manual(values = c('darkgoldenrod1', "#413a6e", 'darkgray')) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5), aspect.ratio = 1) + 
  geom_rect(data = TEs_sub, aes(xmin = V4, xmax = V5, ymin = 0,ymax = 0.2, group = 1), fill = 'gray90', color = 'gray90', alpha = 0.8, linewidth = 0) +
  geom_line(aes(x = window, y = Dxy, color = hap, group = qryGnom)) 
  


# ===== Dxy between J. regia haplotypes =====

DxyhJregSub[hap == 'H', hap := 'G']
DxyhJregSub[hap == 'h', hap := 'g']
plot1 <- ggplot(DxyhJregSub[species == 'Jreg'], aes(x = window, y = Dxy, color = hap, group = qryGnom) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', color = '', title = 'J. regia') +
  geom_rect(data = DxyhJregSub[species == 'Jreg' & qryGnom == 'HJregBNU'], aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line(linewidth = 1) +
  scale_y_continuous(breaks = c(0,0.1,0.2)) +
  
  scale_x_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 3), labels = sprintf("%.3f", seq(31.86, 31.89, length.out = 3))) +
  scale_color_manual(values = c('darkgoldenrod1', "#413a6e", 'darkgray')) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5), 
        aspect.ratio = 0.6, 
        axis.text = element_text(size = 10), 
        legend.position = c(0.9, 0.8),
        plot.title = element_text(face = 'italic', size = 10),
        legend.text = element_text(size = 10, face = 'italic'),
        legend.key.size = unit(0.5, 'cm')) + 
  # NDR1
  annotate("rect", xmin = 31868750,
           xmax = 31869735, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = '#94435b')  +
  #TPPD
  annotate("rect", xmin = 31884475,
           xmax = 31887080, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = 'blue') + 
  # mean
  geom_hline(data = DxyhJregSub[species == 'Jreg' & qryGnom == 'HJregBNU'], 
             aes(yintercept = avgDxy), color = 'darkgray', linetype = 2)  #+
  #geom_segment(aes(x = 31887080, y = 0.195, xend = 31884475, yend = 0.195),
  #             arrow = arrow(length = unit(0.2, "cm")), color = 'blue') + 
  #geom_segment(aes(x = 31868750, y = 0.195, xend = 31869735, yend = 0.195),
  #             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b')



plot1


# ----- plot including outgroup -----


p1_outgrp <- ggplot(DxyhJregSub[species %in% c('Jreg', "Cill")], aes(x = window, y = Dxy, color = hap, group = qryGnom) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', color = '', title = 'J. regia') +
  geom_rect(data = DxyhJregSub[species == 'Jreg' & qryGnom %in% c('HJregBNU', 'CillPaw')], aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line(linewidth = 1) +
  scale_y_continuous(breaks = c(0,0.1,0.2)) +
  
  scale_x_continuous(breaks = seq(31.86e6, 31.89e6, length.out = 3), labels = sprintf("%.3f", seq(31.86, 31.89, length.out = 3))) +
  scale_color_manual(values = c('darkgoldenrod1', "#413a6e", 'turquoise')) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5), 
        aspect.ratio = 0.6, 
        axis.text = element_text(size = 10), 
        legend.position = c(0.9, 1),
        plot.title = element_text(face = 'italic', size = 10),
        legend.text = element_text(size = 10, face = 'italic'),
        legend.key.size = unit(0.5, 'cm')) + 
  # NDR1
  annotate("rect", xmin = 31868750,
           xmax = 31869735, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = '#94435b')  +
  #TPPD
  annotate("rect", xmin = 31884475,
           xmax = 31887080, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = 'blue') + 
  # mean
  geom_hline(data = DxyhJregSub[species == 'Jreg' & qryGnom == 'HJregBNU'], 
             aes(yintercept = avgDxy), color = 'darkgray', linetype = 2)  #+
#geom_segment(aes(x = 31887080, y = 0.195, xend = 31884475, yend = 0.195),
#             arrow = arrow(length = unit(0.2, "cm")), color = 'blue') + 
#geom_segment(aes(x = 31868750, y = 0.195, xend = 31869735, yend = 0.195),
#             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b')

p1_outgrp


# ===== Plot Divergence against hJmanBNU =====
Jman_astart <- 31745000
Jman_aend <- 31790000

DxyhJmanSub <- dxy[refGnom == 'hJmanBNU' &
                     
                     window > Jman_astart &
                     window < Jman_aend ]

plot2 <- ggplot(DxyhJmanSub[qryGnom == 'HJmanNFU' ], aes(x = window, y = Dxy, group = qryGnom) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', title = 'J. mandshurica') +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_hline(aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) +
  
  geom_line(color = "#413a6e", linewidth = 1) +
  
  scale_x_continuous(breaks = seq(31.75e6, 31.79e6, length.out = 3), labels = seq(31.75, 31.79, length.out = 3) ) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5), 
        axis.text = element_text(size = 10),
        plot.title = element_text(face = 'italic', size = 10),
        aspect.ratio = 0.6)  +
  # NDR1
  annotate("rect", xmin = 31760070,
           xmax = 31761382, 
           ymin=0.2, 
           ymax = 0.21, 
           alpha = .5,fill = '#94435b')  +
  #TPPD
  annotate("rect", xmin = 31776425,
           xmax = 31778494, 
           ymin=0.2, 
           ymax = 0.21, 
           alpha = 0.5,fill = 'blue')   
#geom_segment(aes(x = 31760070, y = 0.205, xend = 31761382, yend = 0.205),
#             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
#geom_segment(aes(x = 31778494, y = 0.205, xend = 31776425, yend = 0.205),
#             arrow = arrow(length = unit(0.2, "cm")), color = 'blue') 

plot2


# ------ Plot with pecan as outgroup -----

p2_outgrp <- ggplot(DxyhJmanSub[], aes(x = window, y = Dxy, group = qryGnom, color = hap) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', title = 'J. mandshurica') +
  geom_line(linewidth = 1) +
  
  geom_rect( data = DxyhJmanSub[qryGnom == 'HJmanFNU'], aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_hline( data = DxyhJmanSub[qryGnom == 'HJmanFNU'], aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) +
  
  scale_x_continuous(breaks = seq(31.75e6, 31.79e6, length.out = 3), labels = seq(31.75, 31.79, length.out = 3) ) +
  theme_classic() + 
  scale_color_manual(values = c("#413a6e", 'turquoise')) +
  
  theme(axis.text.x = element_text(vjust = 0.5), 
        axis.text = element_text(size = 10),
        plot.title = element_text(face = 'italic', size = 10),
        aspect.ratio = 0.6,
        legend.position = 'none')  +
  # NDR1
  annotate("rect", xmin = 31760070,
           xmax = 31761382, 
           ymin=0.2, 
           ymax = 0.21, 
           alpha = .5,fill = '#94435b')  +
  #TPPD
  annotate("rect", xmin = 31776425,
           xmax = 31778494, 
           ymin=0.2, 
           ymax = 0.21, 
           alpha = 0.5,fill = 'blue')   
#geom_segment(aes(x = 31760070, y = 0.205, xend = 31761382, yend = 0.205),
#             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
#geom_segment(aes(x = 31778494, y = 0.205, xend = 31776425, yend = 0.205),
#             arrow = arrow(length = unit(0.2, "cm")), color = 'blue') 

p2_outgrp






# ===== Plot Dxy between J. californica haplotypes =====
Jcal_astart <- 30480000
Jcal_aend <- 30536000

DxyhJcalSub <- dxy[refGnom == 'hJcaliAlt' & 
                     
                     window > Jcal_astart & 
                     window < Jcal_aend]


plot3 <- ggplot(DxyhJcalSub[qryGnom == 'HJcaliPrimary'], aes(x = window, y = Dxy, color = hap, group = qryGnom) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', title = 'J. californica') +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line(color = "#413a6e", linewidth = 1) +
  scale_x_continuous(breaks = seq(30.49e6, 30.53e6, length.out = 3), labels = seq(30.49, 30.53, length.out = 3)) +
  scale_y_continuous(breaks = c(0,.1,.2)) +
  
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5), 
        aspect.ratio = 0.6, 
        plot.title = element_text(face = 'italic', size = 10),
        axis.text = element_text(size = 10)) + 
  # NDR1
  annotate("rect", xmin = 30498727,
           xmax = 30500052, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = '#94435b')  +
  # TPP
  annotate("rect", xmin = 30517523,
           xmax = 30520347, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = 'blue')  +

  #geom_segment(aes(x = 30498727, y = 0.195, xend = 30500052, yend = 0.195),
  #             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
  #geom_segment(aes(x = 30520347, y = 0.195, xend = 30517523, yend = 0.195),
  #           arrow = arrow(length = unit(0.2, "cm")), color = 'blue') +
#
  # mean
  geom_hline(data = DxyhJcalSub[qryGnom == 'HJcaliPrimary'], aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) 

plot3


# ------ plot with pecan as outgroup -----
p3_outgrp <- ggplot(DxyhJcalSub, aes(x = window, y = Dxy, color = hap, group = qryGnom) ) + 
  #facet_wrap(~species) +
  #geom_point()  +
  labs(x = "", y = '', title = 'J. californica') +
  geom_rect(data = DxyhJcalSub[qryGnom == 'HJcaliPrimary'], aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(30.49e6, 30.53e6, length.out = 3), labels = seq(30.49, 30.53, length.out = 3)) +
  scale_y_continuous(breaks = c(0,.1,.2)) +
  scale_color_manual(values = c("#413a6e", 'turquoise')) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5), 
        aspect.ratio = 0.6, 
        plot.title = element_text(face = 'italic', size = 10),
        axis.text = element_text(size = 10),
        legend.position = 'none') + 
  # NDR1
  annotate("rect", xmin = 30498727,
           xmax = 30500052, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = '#94435b')  +
  # TPP
  annotate("rect", xmin = 30517523,
           xmax = 30520347, 
           ymin=0.19, 
           ymax = 0.2, 
           alpha = .5,fill = 'blue')  +
  
  #geom_segment(aes(x = 30498727, y = 0.195, xend = 30500052, yend = 0.195),
  #             arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
  #geom_segment(aes(x = 30520347, y = 0.195, xend = 30517523, yend = 0.195),
  #           arrow = arrow(length = unit(0.2, "cm")), color = 'blue') +
  #
  # mean
  geom_hline(data = DxyhJcalSub[qryGnom == 'HJcaliPrimary'], aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) 

p3_outgrp






# ===== Plot in multipanel =====

DxyMultiPanel <- ggarrange(plot1, plot2, plot3, ncol = 1)
DxyMultiPanel


DxyMultiPanel_outgrp <- ggarrange(p1_outgrp, p2_outgrp, p3_outgrp, ncol = 3)
DxyMultiPanel_outgrp  



# ===== Plot divergence against microcarpa =====

DxyHJmicSub <- dxy[refGnom == 'HJmic' & 
                 window > 30.54e6 &  
                 window < 30.58e6]

DxyHJmicSub[species == 'Cycl', species := 'Cpal']
DxyHJmicSub[, species1 := factor(species, levels = c("Jcal", 'Jnig',"Jreg", 'Jsig', 'Jman', 'Pste', 'Cpal', 'Pstr', 'Plon', 'Cill'))]


ggplot(DxyHJmicSub, aes(x = window, y = Dxy)) +
  facet_wrap(~species1) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.95), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_line(aes(color = hap, group = qryGnom), size = 0.8)  +
  scale_color_manual(values = c('darkgoldenrod1', "#413a6e", "darkgray")) +
  #geom_line(data = DxyHJmicSub[hap == 'H'], color = "#413a6e", size = 1) +
  #geom_line(data = DxyHJmicSub[hap == 'h'], color = "darkgoldenrod1", size = 1) +
  labs(x = "H haplotype (Mb)", y = 'Nucleotide divergence', color = '') +
  scale_x_continuous(breaks = seq(30.546e6, 30.58e6, length.out = 3), labels = seq(30.546e6, 30.58e6, length.out = 3)/1e6 ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5),
        legend.position = c(0.9, 0.1), 
        legend.key.size = unit(1, 'cm'), 
        legend.text = element_text(size = 14), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 14)) + 
  #TPPD
  annotate("rect", xmin = 30571569,
           xmax = 30574389, 
           ymin=0, 
           ymax = 0.3, 
           alpha = .3,fill = '#94435b')  +
  #NDR1
  annotate("rect", xmin = 30546511,
           xmax = 30547577, 
           ymin=0, 
           ymax = 0.3, 
           alpha = .3,fill = 'turquoise4')   + 
  #HJ3
  #annotate("rect", xmin = 30.546e6 + 3075,
  #         xmax = 30.546e6 + 5885, 
  #         ymin=0, 
  #         ymax = 0.3, 
  #         alpha = .3,fill = 'tan') + 
  # avg Dxy
  geom_hline(aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) 

 # geom_vline(aes(xintercept = 30.5501e6), color = 'grey', linetype = 2) 
  
  


# ====== Date Estimate =====
JK <- function(p){
  -(3/4)*log(1 - (4/3)*p)
}

# ===== Approach #1. Using within-species alignments and reported substitution rates =====

DxyhJregSub[window > 31.87e6 & window < 31.89e6 & qryGnom == 'HJregBNU']

regL <-  0.092 #(249 bp)
regR <-  0.195 #(133 bp)

DxyhJmanSub[window > 31.76e6 & window < 31.8e6]
manL <- 0.246 #(419 bp)
manR <- 0.254 #(488 bp)

DxyhJcalSub[window > 30.49e6 & window < 30.53e6]
calL <-  0.170 #(446 bp)
calR <- 0.125 # (200 bp)

Hdxy <- (regL + regR + manL + manR + calL + calR)/6

Hdxy_adj <- JK(Hdxy)
# div = 2Tu
# T = div/(2u)

Hdxy_adj / (2*1.5e-9)
Hdxy_adj / (2*2.5e-9)

# ===== Approach 2 =====
JK <- function(p){
  -(3/4)*log(1 - (4/3)*p)
}

# # ---- mandshurica -----
# DxyhJregSub[window > 31.868e6 & window < 31.888e6 & qryGnom == 'HJmanNFU']
# 
# Jman_Chandler_Gloc_div <- 0.24336283
# Jman_Chandler_WG_div <- DxyhJregSub[qryGnom == 'HJmanNFU', avgDxy][1]
# 
# JK(Jman_Chandler_Gloc_div) # 0.2942086
# JK(Jman_Chandler_WG_div) # 0.06529312
# 
# # divergence time of J. man and J. reg 
# # 31.66 (Mu et al. 2020)
# 
# # divergence time of J. man and J. reg 
# # 39.3 (Zhou et al. 2021)
# 
# 
# x1 <- 31.66*(0.06529312/0.2942086)
# x2 <- 39.3*(0.06529312/0.2942086)
# 
# 
# # ----- microcarpa -----
# DxyhJregSub[window > 31.868e6 & window < 31.888e6 & qryGnom == 'HJmic']
# 
# Jmic_Chandler_Gloc_div <- 0.23250564
# Jmic_Chandler_WG_div <- DxyhJregSub[qryGnom == 'HJmic', avgDxy][1]
# 
# JK(Jmic_Chandler_Gloc_div) # 0.2783059
# JK(Jmic_Chandler_WG_div) # 0.06850902
# 
# # divergence time of J. mic and J. reg 
# # 46.71 (Mu et al. 2020)
# 
# # divergence time of J. mic and J. reg 
# # 43.9 (Zhou et al. 2021)
# 
# y1 <- 46.71*(0.2783059/0.06850902)
# y2 <- 43.9*(0.2783059/0.06850902)
# 
# # ----- californica -----
# DxyhJregSub[window > 31.868e6 & window < 31.888e6 & qryGnom == 'HJcaliPrimary']
# 
# Jcal_Chandler_Gloc_div <- 0.27408994
# Jcal_Chandler_WG_div <- DxyhJregSub[qryGnom == 'HJcaliPrimary', avgDxy][1]
# 
# 
# JK(Jcal_Chandler_Gloc_div) # 0.3411332
# JK(Jcal_Chandler_WG_div) # 0.06663864
# 
# 
# z1 <- 46.71*(0.3411332/0.06663864)
# z2 <- 43.9*(0.3411332/0.06663864)
# 
# # divergence time of J. cal and J. reg 
# # 46.71 (Mu et al. 2020)
# # 43.9 (Zhou et al. 2021)


