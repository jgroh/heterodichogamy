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

# -- coverage to BNU assembly
reg2BNU_rawcvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_BNU/", 
                                                  pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))

reg2BNU_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_BNU/", 
                                                      pattern = '*_norm.txt', full.names = T), 
                                           function(x) {
                                             z <- fread(x, col.names = 'avg_cvg')
                                             z[, sample := gsub("_norm.txt", "", basename(x))]
                                             return(z)
                                           }
))
reg2BNU_cvg <- merge(reg2BNU_rawcvg, reg2BNU_cvg_nrm)





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

# ---- calculate coverage in 500bp windows ------
coverage[, window := cut(pos, breaks = seq(31.8e6, 32e6, by = 500), labels = seq(31.8e6, 32e6-500, by = 500), include.lowest =T), by = sample]
coverage[, window := as.numeric(as.character((window)))]
coverage_bins <- coverage[, .(coverage = mean(coverage)), by = .(sample, species, window)]
coverage_bins <- merge(coverage_bins, pheno)


coverage[, window10 := cut(pos, breaks = seq(31.8e6, 32e6, by = 10), labels = seq(31.8e6, 32e6-10, by = 10), include.lowest =T), by = sample]
coverage[, window10 := as.numeric(as.character((window10)))]
coverage_bins10 <- coverage[, .(coverage = mean(coverage)), by = .(sample, species, window10)]
coverage_bins10 <- merge(coverage_bins10, pheno)


#30710000-30817000
winsize <- 1e3
reg2BNU_cvg[, window := cut(position, breaks = seq(30710000, 30817000, by = winsize), labels = seq(30710000, 30817000-winsize, by = winsize), include.lowest =T), by = sample]
reg2BNU_cvg[, window := as.numeric(as.character((window)))]
reg2BNU_cvg_win <- reg2BNU_cvg[, .(coverage = mean(coverage)), by = .(sample, window, avg_cvg)]
reg2BNU_cvg_win <- merge(reg2BNU_cvg_win, pheno[, .(sample, phenotype)])
reg2BNU_cvg_win[, nrm_cvg := coverage/avg_cvg]

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
#31883560-31884214
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
coverage_bins[phenotype == 'protandrous', genotype := 'gg']
coverage_bins[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'Gg']
coverage_bins[sample == 'JG0026', genotype := 'GG']

coverage_bins[, cvg_nrm := coverage/mean(coverage), by = sample]


coverage_bins10[phenotype == 'protandrous', genotype := 'gg']
coverage_bins10[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'Gg']
coverage_bins10[sample == 'JG0026', genotype := 'GG']
coverage_bins10[, cvg_nrm := coverage/mean(coverage), by = sample]

coverage1 <- merge(coverage, pheno)[species == 'regia']
coverage1[phenotype == 'protandrous', genotype := 'gg']
coverage1[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'Gg']
coverage1[sample == 'JG0026', genotype := 'GG']
coverage1[, cvg_nrm := coverage/mean(coverage), by = sample]

plot2 <- ggplot(coverage_bins[species == 'regia' & window >= 31.865e6 & window <= 31.895e6],
       aes(x = window, y = cvg_nrm, group = sample, color = genotype)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("tan", "maroon", 'turquoise4')) +
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
  labs(x = 'Position (Mb)', y =  'Normalized\nread depth', color = '') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.9, 0.8),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 12, face = 'italic')
  ) + 
  geom_segment(aes(x = 31883277, xend =  31884479, y = 4, yend = 4))
plot2

plot2_no_leg <- ggplot(coverage_bins[species == 'regia' & window >= 31.865e6 & window <= 31.895e6],
                aes(x = window, y = cvg_nrm, color = genotype, group = sample)) +
  geom_line(linewidth = 0.8) + 
  scale_color_manual(values = c("tan", "turquoise4", 'maroon')) +
  theme_classic() + 
 # scale_x_continuous(n.breaks = 10, labels = NULL) +
  scale_x_continuous(breaks = c(31.87e6, 31.88e6, 31.89e6), labels = c(31.87, 31.88, 31.89)) +
  
  labs(x = '', y =  'Normalized read depth', color = '') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.09, 0.9),
        
        plot.title = element_text(vjust = -10, hjust = 0.1),
        #legend.position = 'none',
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 12)
        
  )

plot2_no_leg



# ----- coverage to Chandler, finer scale --------

# Chandler TEs
cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
hJregTE <- fread("~/workspace/heterodichogamy/Juglans_genome_assemblies/Jregia/Chandler/JregiaV2.fa.mod.EDTA.TEanno.gff3", skip = "NC_", col.names = cnames)

ggplot(coverage1[species == 'regia' & pos >= 31.883e6 & pos <= 31.886e6],
       #species == 'regia' & pos >= 31.88447e6 & pos <= 31.88449e6], # showing boundary in detail
                       ) +
  #geom_point() +
  geom_line(aes(x = pos, y = cvg_nrm, color = genotype, group = sample), linewidth = 0.8) + 
  scale_color_manual(values = c("tan", "turquoise4", 'maroon')) +
  theme_classic() + 
  #scale_x_continuous(breaks = c(31.87e6, 31.88e6, 31.89e6), labels = c(31.87, 31.88, 31.89)) +
  
  labs(x = '', y =  'Normalized read depth', color = '') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.09, 0.9),
        
        plot.title = element_text(vjust = -10, hjust = 0.1),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 12)
  ) + 
  #geom_rect(data = hJregTE[Chr == 'NC_049911.1' & Start > 31.883e6 & End < 31.885e6 ],
  #          aes(ymin = -Inf, ymax = Inf, xmin = Start, xmax = End), fill = 'lightgray') + 
  # UTR
  annotate("rect", xmin = 31884269, xmax = 31885000, ymin = 4, ymax = 6,
         alpha = .1,fill = '#94435b') 
  # starting from last CDS
  #annotate("rect", xmin = 31884478, xmax = 31884490, ymin = 4, ymax = 6,
  #         alpha = .1,fill = '#94435b') 

# full TPP coordinates
#31884269        31887072







# ----- coverage to BNU----- 
reg2BNU_cvg_win[phenotype == 'protandrous', genotype := 'gg']
reg2BNU_cvg_win[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'Gg']
reg2BNU_cvg_win[sample == 'JG0026', genotype := 'GG']
BNUst <- 30746000
BNUen <- 30805000
  
cvg_plt2BNU <- ggplot(reg2BNU_cvg_win[genotype != '??' & window > BNUst & window < BNUen],
                  aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  geom_line(data = reg2BNU_cvg_win[ genotype != '??' & window > BNUst & window < BNUen], linewidth = 0.8, alpha = 0.9)  +
  
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  scale_color_manual( values = c('tan', 'maroon', 'turquoise4')) +
  scale_x_continuous(limits = c(BNUst, BNUen), breaks = seq(30.75e6, 30.79e6, length.out = 3), labels = seq(30.75, 30.79, length.out = 3)) +
  
  theme_classic() + 
  #scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized\nread depth', color = '', title = '') +
  theme(
    #aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.position = c(0.9, 1.1),
    legend.text = element_text(size = 10, face = 'italic'),
    legend.key.size = unit(.4, 'cm')
  ) +
annotate("rect", xmin = 30759890, xmax = 30760855, ymin = 1.8, ymax = 2,
         alpha = .1,fill = '#94435b') + 
  geom_segment(aes(x = 30759890, y = 1.9, xend = 30760855, yend = 1.9),
               arrow = arrow(length = unit(0.2, "cm")), color = '#94435b') +
  annotate("rect", xmin = 30784032, xmax = 30786419, ymin = 1.8, ymax = 2,
           alpha = .1,fill = 'blue') + 
  geom_segment(aes(x = 30786419, y = 1.9, xend = 30784032, yend = 1.9),
               arrow = arrow(length = unit(0.2, "cm")), color = 'blue') #+
cvg_plt2BNU

ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/Jregia_cvg2BNU.pdf', plot = cvg_plt2BNU, width = 7, height = 2.5, units = 'in')











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
