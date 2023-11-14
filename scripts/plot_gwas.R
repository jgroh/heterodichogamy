library(data.table)
library(ggplot2)
library(cowplot)

# --- read P vals
rad <- fread("~/workspace/heterodichogamy/data/heterodichogamy_tassel_MLM_heterodichogamy_stats_trim.txt")
wgs <- fread("~/workspace/heterodichogamy/regia/output/Jregia_founders_gemma.assoc.txt")

wgs[, N:= seq_len(.N)]

# --- read coverage
regia_coverage <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_chr11_31.8-32Mb/", full.names = T),
       function(x){
         z <- fread(x, select = 2:3, col.names = c("pos", "coverage"))
         z[, sample := gsub(".txt.gz", "", basename(x))]
         return(z)
       }))
regia_coverage[, species := 'regia']


hindsii_coverage <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/hindsii/results/coverage_chr11_31.8-32Mb/", full.names = T),
                                   function(x){
                                     z <- fread(x, select = 2:3, col.names = c("pos", "coverage"))
                                     z[, sample := gsub(".txt.gz", "", basename(x))]
                                     return(z)
                                   }))
hindsii_coverage[, species := 'hindsii']


microcarpa_coverage <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/microcarpa/results/coverage_chr11_31.8-32Mb/", full.names = T),
                                     function(x){
                                       z <- fread(x, select = 2:3, col.names = c("pos", "coverage"))
                                       z[, sample := gsub(".txt.gz", "", basename(x))]
                                       return(z)
                                     }))
microcarpa_coverage[, species := 'microcarpa']

nigra_coverage <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/nigra/results/coverage_chr11_31.8-32Mb/", full.names = T),
                                        function(x){
                                          z <- fread(x, select = 2:3, col.names = c("pos", "coverage"))
                                          z[, sample := gsub(".txt.gz", "", basename(x))]
                                          return(z)
                                        }))
nigra_coverage[, species := 'nigra']

coverage <- rbindlist(list(regia_coverage, hindsii_coverage, microcarpa_coverage, nigra_coverage))



# --- read phenotypes 
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("sample", 'phenotype', 'species'))
pheno[sample == 'JG0026', phenotype := 'protogynous (homozygous)']


# --- read dxy
dxy <- fread("~/workspace/heterodichogamy/regia/results/pixy/protandrous_vs_protogynous_dxy.txt")


# ----- read gene tracks
genes <- fread("~/workspace/heterodichogamy/H_locus_structure/Chandler_genes_31.865-31.895Mb.BED", 
               select = c(2,3,6), col.names = c("start", "end", 'strand'))

# ---- calculate coverage in 500bp windows
coverage[, window := cut(pos, breaks = seq(31.8e6, 32e6, by = 500), labels = seq(31.8e6, 32e6-500, by = 500), include.lowest =T), by = sample]
coverage[, window := as.numeric(as.character((window)))]
coverage_bins <- coverage[, .(coverage = mean(coverage)), by = .(sample, species, window)]
coverage_bins <- merge(coverage_bins, pheno)


# ===== make plots ===== 
 
# ---- gwas -----

# whole genome (takes a minute to plot)
# *note to fix x axis for final plot
# for the graphic to load in reasonable time subset data except in region of interest
plt_GWAS_data <- rbind(wgs[chr != 'NC_049911.1' | (chr == 'NC_049911.1' & ps < 31e6) | (chr == 'NC_049911.1' & ps > 32e6)][seq(1, .N, by = 10)], wgs[chr == 'NC_049911.1' & ps > 31e6 & ps < 32e6])

ggplot(plt_GWAS_data, aes(x = N, y = -log10(p_wald), color = chr)) + 
         geom_point(size = 1) + 
         scale_color_manual(values = rep(c("black", "grey"), 8)) + 
  theme_classic() + 
  labs(y = expression(-log[10](P)), 
       x = '') +
  theme(aspect.ratio = .3,
        plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2),
        #legend.position = c(0.15, 0.85), 
        legend.position = 'none')

# save as 7.5 x 5.5 

# just chromosome 11
#ggplot(wgs[CHR=='NC_049911.1'], aes(x = BP, y = -log10(P))) + 
#  geom_point() 

# finer scale
plot1 <- ggplot(wgs[chr=='NC_049911.1' & ps > 3.1865e7 & ps < 3.1895e7], aes(x = ps, y = -log10(p_wald))) + 
  geom_point() + 
  theme_classic() + 
  labs(x = '', y = expression(-log[10] (P))) +
  scale_x_continuous(breaks = c(31.87e6, 31.88e6, 31.89e6), labels = NULL)  +
  theme(aspect.ratio = 0.25,
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  ) + 
  annotate("rect", xmin = genes$start[1], xmax = genes$end[1], ymin = 16, ymax = 18,
           alpha = .1,fill = '#94435b') + 
  geom_segment(aes(x = genes$start[1], y = 17, xend = genes$end[1], yend = 17),
               arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
  annotate("rect", xmin = genes$start[2], xmax = genes$end[2], ymin = 16, ymax = 18,
           alpha = .1,fill = 'blue') + 
  geom_segment(aes(x = genes$end[2], y = 17, xend = genes$start[2], yend = 17),
               arrow = arrow(length = unit(0.2, "cm")), color = 'blue') #+
  #annotate('text', x = (genes$end[4] + genes$start[4])/2, y = 14, label='***') 
  
  
  
plot1

# ----- coverage -----
coverage_bins[phenotype == 'protandrous', genotype := 'hh']
coverage_bins[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'Hh']
coverage_bins[sample == 'JG0026', genotype := 'HH']

plot2 <- ggplot(coverage_bins[species == 'regia' & window >= 31.865e6 & window <= 31.895e6],
       aes(x = window, y = coverage, group = sample, color = genotype)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tan", "turquoise4", 'maroon')) +
  #geom_line(data = coverage_bins[species == 'regia' & genotype == 'hh' & window >= 31.865e6 & window <= 31.895e6], 
  #          aes(x = window, y = coverage), 
  #          alpha = 0.9, linewidth = 0.8, color = 'tan') + 
  #geom_line(data = coverage_bins[species == 'regia' & genotype == 'Hh' & window >= 31.865e6 & window <= 31.895e6], 
  #          aes(x = window, y = coverage),
  #          alpha = 0.9, linewidth = 0.8, color = 'turquoise4') + 
  #geom_line(data = coverage_bins[species == 'regia' & genotype == 'HH' & window >= 31.865e6 & window <= 31.895e6], 
  #          aes(x = window, y = coverage),
  #          alpha = 0.9, linewidth = 0.8, color = 'maroon') + 
  
  #scale_color_manual(values = c("#cfb582","#28bbd1",'#94435b')) + 
  theme_classic() + 
  scale_x_continuous(breaks = c(31.87e6, 31.88e6, 31.89e6), labels = c(31.87, 31.88, 31.89)) +
  labs(x = 'Position (Mb)', y =  'Read depth', color = '') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.1, 0.9),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 12)
  )
plot2

plot2_no_leg <- ggplot(coverage_bins[species == 'regia' & window >= 31.86e6 & window <= 31.9e6],
                aes(x = window, y = coverage, color = phenotype, group = sample)) +
  geom_line(linewidth = 1) + 
  scale_color_manual(values = c("#cfb582","#28bbd1",'#94435b')) + 
  theme_classic() + 
  scale_x_continuous(n.breaks = 10, labels = NULL) +
  labs(x = '', y =  'Read depth', color = '', title = 'regia') +
  theme(aspect.ratio = 0.25,
        plot.title = element_text(vjust = -10, hjust = 0.1),
        legend.position = 'none',
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  )

# ----- Dxy -----

plot3 <- ggplot(dxy[window_pos_1> 31.86e6 & window_pos_2 < 31.905e6 ], 
       aes(x = window_pos_1, y = avg_dxy)) + geom_line(linewidth = 1) + 
  scale_x_continuous(breaks = seq(31.86e6, 31.9e6, length.out = 5), labels = seq(31.86e6, 31.9e6, length.out = 5)/1e6 ) +
  theme_classic() + 
  labs(x = 'Position on chromosome 11 (Mb)', y = expression(D[XY]) , color = '') + 
  theme(aspect.ratio = 0.25,
        plot.margin = unit(c(1,1,1,1), 'lines'),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  )

plot3
plot_grid(plot1, plot2, plot3, ncol=1, align = 'hv')


hindsii_coverage <- ggplot(coverage_bins[species == 'hindsii' & window >= 31.86e6 & window <= 31.9e6],
aes(x = window, y = coverage, group = sample, color = phenotype)) +
  #scale_color_manual(values = 'gray') +
  geom_line(linewidth = 1) + 
  #scale_color_manual(values = c("#cfb582","#28bbd1",'#94435b')) + 
  theme_classic() + 
  scale_x_continuous(n.breaks = 10, labels = NULL) +
  labs(x = '', y =  'Read depth', title = 'hindsii') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.9,.8),
       # plot.margin = unit(c(0,0,0,0), "lines"),
        plot.title = element_text(vjust = -10, hjust = 0.1, size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  )

hindsii_coverage
nigra_coverage <- ggplot(coverage_bins[species == 'nigra' & window >= 31.86e6 & window <= 31.9e6],
                           aes(x = window, y = coverage, color = phenotype, group = sample)) +
  scale_color_manual(values = c('darkred', 'gray')) +
  geom_line(linewidth = 1) + 
  #scale_color_manual(values = c("#cfb582","#28bbd1",'#94435b')) + 
  theme_classic() + 
  scale_x_continuous(n.breaks = 10, labels = NULL) +
  labs(x = '', y =  'Read depth', title = 'nigra') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.9,0.8),
        #plot.margin = unit(c(0,0,0,0), "lines"),
        plot.title = element_text(vjust = -10, hjust = 0.1, size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  )
nigra_coverage

microcarpa_coverage <- ggplot(coverage_bins[species == 'microcarpa' & window >= 31.86e6 & window <= 31.9e6],
                           aes(x = window, y = coverage, color = phenotype, group = sample)) +
  geom_line(linewidth = 1) + 
  scale_color_manual(values = c('darkred', 'gray')) +
  theme_classic() + 
  scale_x_continuous(n.breaks = 10, labels = seq(31.86e6, 31.9e6, length.out=9)/1e6 ) +
  labs(x = 'Position on chromosome 11 (Mb)', y =  'Read depth', title = 'microcarpa') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.9,0.8),
        #plot.margin = unit(c(-1,0,0,0), "lines"),
        plot.title = element_text(vjust = -10, hjust = 0.1,size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
  )

microcarpa_coverage

plot(x = coverage_bins[species == 'microcarpa', coverage], y = coverage_bins[species == 'nigra', coverage])

# ----- compare RAD and WGS GWAS peaks

chr11_rad <- rad[Chr == 11 & Trait == 'protog']
chr11_rad[, data := 'RAD']
chr11_wgs <- wgs[CHR == 'NC_049911.1']
chr11_wgs[, data := 'WGS']

# combine

all_P <- rbind(chr11_rad[, .(Pos, P = p, data)],
      chr11_wgs[, .(Pos = BP, P, data)])

ggplot(all_P, aes(x = Pos, y = -log10(P), color = data)) + 
  geom_point()  + 
  scale_color_manual(values = c("gray", "black")) + 
  labs(x = "Position Chr 11", color = '') +
  theme_classic() + 
  theme(legend.position = c(0.1,0.9))

ggplot(all_P[Pos > 3.15e7 & Pos < 3.2e7 & -log10(P) > 2], aes(x = Pos, y = -log10(P), color = data)) + 
  geom_point()  + 
  scale_color_manual(values = c("gray", "black"))

ggplot(all_P[Pos > 3.18e7 & Pos < 3.2e7], aes(x = Pos, y = -log10(P), color = data)) + 
  geom_point()  + 
  scale_color_manual(values = c("red", "black")) + 
  theme_classic()

ggplot(all_P[Pos > 3.188e7 & Pos < 3.189e7], aes(x = Pos, y = -log10(P))) + 
  geom_point()  

all_P[data == 'WGS'][P == min(P)]
