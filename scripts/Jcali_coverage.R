library(data.table)

library(data.table)
library(ggplot2)
library(cowplot)


Jcal_pstart <- 31348000
Jcal_pend <- 31400000

Jcal_astart <- 30490000
Jcal_aend <- 30536000

pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", col.names = c("sample", 'phenotype', 'species'))
pheno[phenotype == 'protogynous' & sample != 'JG0026', genotype := 'H?']
pheno[sample == 'JG0026', genotype := 'HH']
pheno[phenotype == 'protandrous', genotype := 'hh']
pheno[phenotype == 'unknown', genotype := '??']




# ===== Read coverage =====

# ----- hindsii
# -- to cali primary assembly
hindsii2cali_p_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/hindsii/results/coverage_Jcali_primary//", 
                                                  pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
hindsii2cali_p_cvg[, species := 'hindsii']
hindsii2cali_p_cvg[, reference := 'primary']

hindsii2cali_p_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/hindsii/results/coverage_Jcali_primary/", 
                  pattern = '*_norm.txt', full.names = T), 
       function(x) {
         z <- fread(x, col.names = 'avg_cvg')
         z[, sample := gsub("_norm.txt", "", basename(x))]
         return(z)
       }
))
hindsii2cali_p_cvg_nrm[, species := 'hindsii']
hindsii2cali_p_cvg_nrm[, reference := 'primary']

# -- to cali alternate assembly
hindsii2cali_alt_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/hindsii/results/coverage_Jcali_alt//", 
                                                    pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))
hindsii2cali_alt_cvg[, species := 'hindsii']
hindsii2cali_alt_cvg[, reference := 'alternate']

hindsii2cali_alt_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/hindsii/results/coverage_Jcali_alt/", 
                                                      pattern = '*_norm.txt', full.names = T), 
                                           function(x) {
                                             z <- fread(x, col.names = 'avg_cvg')
                                             z[, sample := gsub("_norm.txt", "", basename(x))]
                                             return(z)
                                           }
))
hindsii2cali_alt_cvg_nrm[, species := 'hindsii']
hindsii2cali_alt_cvg_nrm[, reference := 'alternate']


# ----- microcarpa
# -- to cali primary assembly
mic2cali_p_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/microcarpa/results/coverage_Jcali_primary//", 
                                              pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))
mic2cali_p_cvg[, species := 'microcarpa']
mic2cali_p_cvg[, reference := 'primary']

mic2cali_p_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/microcarpa/results/coverage_Jcali_primary/", 
                                                        pattern = '*_norm.txt', full.names = T), 
                                             function(x) {
                                               z <- fread(x, col.names = 'avg_cvg')
                                               z[, sample := gsub("_norm.txt", "", basename(x))]
                                               return(z)
                                             }
))
mic2cali_p_cvg_nrm[, species := 'microcarpa']
mic2cali_p_cvg_nrm[, reference := 'primary']

# -- to cali alternate assembly
mic2cali_alt_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/microcarpa/results/coverage_Jcali_alt//", 
                                                pattern = '*.txt.gz', full.names = T),
                                   function(x){
                                     z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                     z[, sample := gsub(".txt.gz", "", basename(x))]
                                     return(z)
                                   }))
mic2cali_alt_cvg[, species := 'microcarpa']
mic2cali_alt_cvg[, reference := 'alternate']

mic2cali_alt_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/microcarpa/results/coverage_Jcali_alt/", 
                                                    pattern = '*_norm.txt', full.names = T), 
                                         function(x) {
                                           z <- fread(x, col.names = 'avg_cvg')
                                           z[, sample := gsub("_norm.txt", "", basename(x))]
                                           return(z)
                                         }
))
mic2cali_alt_cvg_nrm[, species := 'microcarpa']
mic2cali_alt_cvg_nrm[, reference := 'alternate']

# ----- nigra
# -- to cali primary assembly
nigra2cali_p_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/nigra/results/coverage_Jcali_primary//", 
                                                pattern = '*.txt.gz', full.names = T),
                                   function(x){
                                     z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                     z[, sample := gsub(".txt.gz", "", basename(x))]
                                     return(z)
                                   }))
nigra2cali_p_cvg[, species := 'nigra']
nigra2cali_p_cvg[, reference := 'primary']

nigra2cali_p_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/nigra/results/coverage_Jcali_primary/", 
                                                    pattern = '*_norm.txt', full.names = T), 
                                         function(x) {
                                           z <- fread(x, col.names = 'avg_cvg')
                                           z[, sample := gsub("_norm.txt", "", basename(x))]
                                           return(z)
                                         }
))

nigra2cali_p_cvg_nrm[, species := 'nigra']
nigra2cali_p_cvg_nrm[, reference := 'primary']


# -- to cali alternate assembly
nigra2cali_alt_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/nigra/results/coverage_Jcali_alt//", 
                                                  pattern = '*.txt.gz', full.names = T),
                                     function(x){
                                       z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                       z[, sample := gsub(".txt.gz", "", basename(x))]
                                       return(z)
                                     }))
nigra2cali_alt_cvg[, species := 'nigra']
nigra2cali_alt_cvg[, reference := 'alternate']

nigra2cali_alt_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/nigra/results/coverage_Jcali_alt/", 
                                                    pattern = '*_norm.txt', full.names = T), 
                                         function(x) {
                                           z <- fread(x, col.names = 'avg_cvg')
                                           z[, sample := gsub("_norm.txt", "", basename(x))]
                                           return(z)
                                         }
))
nigra2cali_alt_cvg_nrm[, species := 'nigra']
nigra2cali_alt_cvg_nrm[, reference := 'alternate']


# ----- regia
# -- to cali primary assembly
reg2cali_p_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Jcali_primary//", 
                                              pattern = '*.txt.gz', full.names = T),
                                   function(x){
                                     z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                     z[, sample := gsub(".txt.gz", "", basename(x))]
                                     return(z)
                                   }))
reg2cali_p_cvg[, species := 'regia']
reg2cali_p_cvg[, reference := 'primary']

reg2cali_p_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Jcali_primary/", 
                                                    pattern = '*_norm.txt', full.names = T), 
                                         function(x) {
                                           z <- fread(x, col.names = 'avg_cvg')
                                           z[, sample := gsub("_norm.txt", "", basename(x))]
                                           return(z)
                                         }
))
reg2cali_p_cvg_nrm[, species := 'regia']
reg2cali_p_cvg_nrm[, reference := 'primary']

# -- to cali alternate assembly
reg2cali_alt_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Jcali_alt//", 
                                                pattern = '*.txt.gz', full.names = T),
                                     function(x){
                                       z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                       z[, sample := gsub(".txt.gz", "", basename(x))]
                                       return(z)
                                     }))
reg2cali_alt_cvg[, species := 'regia']
reg2cali_alt_cvg[, reference := 'alternate']

reg2cali_alt_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Jcali_alt/", 
                                                  pattern = '*_norm.txt', full.names = T), 
                                       function(x) {
                                         z <- fread(x, col.names = 'avg_cvg')
                                         z[, sample := gsub("_norm.txt", "", basename(x))]
                                         return(z)
                                       }
))
reg2cali_alt_cvg_nrm[, species := 'regia']
reg2cali_alt_cvg_nrm[, reference := 'alternate']



all_cali_cvg <- rbindlist(list(reg2cali_p_cvg, reg2cali_alt_cvg, mic2cali_p_cvg, mic2cali_alt_cvg, nigra2cali_p_cvg, nigra2cali_alt_cvg, hindsii2cali_p_cvg, hindsii2cali_alt_cvg))
all_cali_cvg_nrm <- rbindlist(list(reg2cali_p_cvg_nrm, reg2cali_alt_cvg_nrm, mic2cali_p_cvg_nrm, mic2cali_alt_cvg_nrm, nigra2cali_p_cvg_nrm, nigra2cali_alt_cvg_nrm, hindsii2cali_p_cvg_nrm, hindsii2cali_alt_cvg_nrm))

all_cali_cvg <- merge(all_cali_cvg, all_cali_cvg_nrm)

# ---- calculate coverage in 500bp windows
all_cali_cvg[reference == 'primary', window := cut(position, breaks = seq(Jcal_pstart, Jcal_pend+500, by = 500), labels = seq(Jcal_pstart, Jcal_pend, by = 500), include.lowest =T), by = sample]
all_cali_cvg[reference == 'alternate', window := cut(position, breaks = seq(Jcal_astart, Jcal_aend+500, by = 500), labels = seq(Jcal_astart, Jcal_aend, by = 500), include.lowest =T), by = sample]

cvg500bp <- all_cali_cvg[, .(coverage = mean(coverage)), by = .(sample, species, window, reference, avg_cvg)]
cvg500bp[, window := as.numeric(as.character((window)))]

cvg500bp <- merge(cvg500bp, pheno)
cvg500bp[, nrm_cvg := coverage/avg_cvg]


# one nigra sample looks HH
cvg500bp[species == 'nigra' & nrm_cvg > 4]
cvg500bp <- cvg500bp[sample != 'JHIN07'] # sample appears to be a duplicate


cvg500bp[, genotype := factor(genotype, levels = c("??", "hh", "H?", 'HH'))]
cvg500bp[, species := factor(species, levels = c("regia", "microcarpa", "hindsii", "nigra"))]
cvg500bp[, reference := factor(reference, levels = c("primary", "alternate"))]



# ---- calculate coverage in 50bp windows
all_cali_cvg[reference == 'primary', window50 := cut(position, breaks = seq(Jcal_pstart, Jcal_pend+50, by = 50), labels = seq(Jcal_pstart, Jcal_pend, by = 50), include.lowest =T), by = sample]
all_cali_cvg[reference == 'alternate', window50 := cut(position, breaks = seq(Jcal_astart, Jcal_aend+50, by = 50), labels = seq(Jcal_astart, Jcal_aend, by = 50), include.lowest =T), by = sample]

cvg50bp <- all_cali_cvg[, .(coverage = mean(coverage)), by = .(sample, species, window50, reference, avg_cvg)]
cvg50bp[, window50 := as.numeric(as.character((window50)))]

cvg50bp <- merge(cvg50bp, pheno)
cvg50bp[, nrm_cvg := coverage/avg_cvg]


cvg50bp[species == 'nigra' & nrm_cvg > 4]
cvg50bp <- cvg50bp[sample != 'JHIN07'] # sample appears to be a duplicate

cvg50bp[, genotype := factor(genotype, levels = c("??", "hh", "H?", 'HH'))]
cvg50bp[, species := factor(species, levels = c("regia", "microcarpa", "hindsii", "nigra"))]
cvg50bp[, reference := factor(reference, levels = c("primary", "alternate"))]




# ==== Plot coverage 500 bp windows =====

cvg500bp
cvg500bp[reference == "primary", refHap := 'H assembly']
cvg500bp[reference == 'alternate', refHap := 'h assembly']
cvg500bp[reference == 'alternate', TPPstart := 30517523]
cvg500bp[reference == 'alternate', TPPend := 30520373]
cvg500bp[reference == 'primary', TPPstart := 31380787]
cvg500bp[reference == 'primary', TPPend := 31383637]


#J. cali alt TPP coords 30517523-30520373
#J. cali primary TPP coords 31380787-31383637
ggplot(cvg500bp[genotype != '??'],
                aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  facet_grid(species~refHap, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  geom_line(data = cvg500bp[ genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  scale_color_manual( values = c('tan', 'turquoise4', 'maroon')) +
  
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized depth', color = '', title = '') +
  theme(aspect.ratio = 0.25,
        #legend.position = c(.49,0.93),
        legend.key.size = unit(0.7, "cm"),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(0, 10, 0, 0),
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
 #       axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
  #      axis.title.x = element_text(size = 14, vjust = -2),
  #      legend.text = element_text(size = 10),
  ) + 
  geom_rect(aes(xmin = TPPend, xmax = TPPstart, ymin=3, ymax = 4), color = 'NA', fill = 'lightgray') +
  geom_segment(aes(x = TPPend, y = 3.5, xend = TPPstart, yend = 3.5),
                    arrow = arrow(length = unit(0.15, "cm")), color = 'black', linewidth = 0.5) 
#annotate("rect", xmin = 30517523, xmax = 30520373, ymin = 3, ymax = 4,
          # alpha = .1,fill = '#94435b')  

cvg500bp[species == 'hindsii' & reference == 'alternate' & window > 30515000 & window < 30520000 & genotype == 'hh' ]

cvg500bp[species == 'hindsii' & reference == 'primary' & window > 31370000 & window < 31372500 & genotype == 'hh' ]







#### ascertain likely protandrous individuals for targeted short read assembly
ggplot(cvg500bp[reference == 'primary' & species %in% c("microcarpa", 'hindsii') & window > 31365000 & window < 31380000],
       aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  facet_grid(species~reference, scales = 'free_x') + 
  geom_line(data = cvg500bp[reference == 'primary' & species %in% c('microcarpa', 'hindsii') & genotype == '??' &  window > 31365000 & window < 31380000], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized depth', color = '', title = '') +
  theme(
    #legend.position = c(.49,0.93),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(0, 10, 0, 0),
    axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
    #       axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12),
    #      axis.title.x = element_text(size = 14, vjust = -2),
    #      legend.text = element_text(size = 10),
  )

ids <- cvg500bp[reference == 'primary' & 
           species %in% c("microcarpa", 'hindsii') & 
           window > 31370000 & window < 31375000 &
           nrm_cvg < 0.1, unique(unlist(sample))]

cvg500bp[ sample %in% ids, suspected_genotype := 'hh']
#### now re plot 
ggplot(cvg500bp[reference == 'primary' & species %in% c("microcarpa", 'hindsii') & window > 31380000 & window < 31390000],
       aes(x = window, y = coverage, group = sample, color = suspected_genotype)) +
  facet_grid(species~reference, scales = 'free_x') + 
  geom_line(data = cvg500bp[reference == 'primary' & species %in% c('microcarpa', 'hindsii') & genotype == '??' & window > 31380000 & window < 31390000], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized depth', color = '', title = '') +
  theme(
    #legend.position = c(.49,0.93),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(0, 10, 0, 0),
    axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
    #       axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12),
    #      axis.title.x = element_text(size = 14, vjust = -2),
    #      legend.text = element_text(size = 10),
  )




# now replot








cvg500bp[species == 'regia' & genotype != 'HH' & 
           genotype != 'hh', genotype := "Hh"]







ggplot(cvg500bp[reference == 'alternate' & species == 'regia'],
       aes(x = window, y = nrm_cvg, group = sample)) +
  #facet_grid(species~reference, scales = 'free_x') + 
  geom_line(data = cvg500bp[reference == 'alternate' & species == 'regia' & genotype == 'hh'], 
            alpha = 0.9, col = "#cfb582", linewidth = 1)  +
  geom_line(data = cvg500bp[reference == 'alternate' & species == 'regia' & genotype == 'Hh'], 
            alpha = 0.9, col = "#28bbd1", linewidth = 1)  +
  geom_line(data = cvg500bp[reference == 'alternate' & species == 'regia' & genotype == 'HH'], 
            alpha = 0.9, col = '#94435b', linewidth = 1)  +
  theme_classic() + 
  labs(x = 'Position on J. cali alt chromosome 11 (Mb)', y =  'Normalized depth', color = '') + 
  scale_x_continuous(breaks = seq(30.49e6, 30.53e6, length.out = 5), labels = seq(30.49e6, 30.53e6, length.out = 5)/1e6 ) +
  theme(aspect.ratio = 0.25,
        legend.position = c(.2,0.9),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 10),
  )

ggplot(cvg500bp[reference == 'primary' & species == 'regia'],
       aes(x = window, y = nrm_cvg, group = sample)) +
  #facet_grid(species~reference, scales = 'free_x') + 
  geom_line(data = cvg500bp[reference == 'primary' & species == 'regia' & genotype == 'hh'], 
            alpha = 0.9, col = "#cfb582", linewidth = 1)  +
  geom_line(data = cvg500bp[reference == 'primary' & species == 'regia' & genotype == 'Hh'], 
            alpha = 0.9, col = "#28bbd1", linewidth = 1)  +
  geom_line(data = cvg500bp[reference == 'primary' & species == 'regia' & genotype == 'HH'], 
            alpha = 0.9, col = '#94435b', linewidth = 1)  +
  theme_classic() + 
  labs(x = 'Position on J. cali v1 chromosome 11 (Mb)', y =  'Normalized depth', color = '') + 
  #scale_x_continuous(breaks = seq(30.49e6, 30.53e6, length.out = 5), labels = seq(30.49e6, 30.53e6, length.out = 5)/1e6 ) +
  scale_x_continuous(breaks = seq(31.35e6, 31.4e6, length.out = 6), labels = seq(31.35e6, 31.4e6, length.out = 6)/1e6 ) +
  
  theme(aspect.ratio = 0.25,
        legend.position = c(.2,0.9),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 10),
  )












# which nigra individual is homozygote? Hay
# cvg500bp[species == 'nigra' & reference == 'alternate' & nrm_cvg == 0]



# ===== Plot coverage in 50 bp windows
#31368050-31379890
cvg50bp[, window50 := as.numeric(window50)]
ggplot(cvg50bp[reference == 'primary' & window50 > 31368050 & window50 < 31379890],
       aes(x = window50, y = nrm_cvg, group = sample, color = genotype)) +
  facet_wrap(~species, scales = 'free_x') + 
  geom_line(data = cvg50bp[reference == 'primary' & genotype == '??' & window50 > 31368050 & window50 < 31379890], linewidth = 0.5, alpha = 0.9)  +
  geom_line(data = cvg50bp[reference == 'primary' & genotype != '??' & window50 > 31368050 & window50 < 31379890], linewidth = 0.5, alpha = 0.9)  +
  scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized depth', color = '', title = '') +
  theme(
    #legend.position = c(.49,0.93),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(0, 10, 0, 0),
    axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
    #       axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12),
    #      axis.title.x = element_text(size = 14, vjust = -2),
    #      legend.text = element_text(size = 10),
  )



