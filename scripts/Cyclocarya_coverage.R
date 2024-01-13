library(data.table)
library(ggplot2)

samples <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/wgs_samples.tsv")
setnames(samples, 'run', 'sample')

# -- coverage to 2PA2023
cyclo_2PA_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Cyclo2PA//", 
                                                   pattern = '*.txt.gz', full.names = T),
                                        function(x){
                                          z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                          z[, sample := gsub(".txt.gz", "", basename(x))]
                                          return(z)
                                        }))


cyclo_2PA_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Cyclo2PA/", 
                                                   pattern = '*_norm.txt', full.names = T), 
                                        function(x) {
                                          z <- fread(x, col.names = 'avg_cvg')
                                          z[, sample := gsub("_norm.txt", "", basename(x))]
                                          return(z)
                                        }
))
cyclo_2PA_cvg <- merge(cyclo_2PA_raw_cvg, cyclo_2PA_cvg_nrm)




# -- coverage to 2PG2023
cyclo_2PG_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Cyclo2PG//", 
                                                 pattern = '*.txt.gz', full.names = T),
                                      function(x){
                                        z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                        z[, sample := gsub(".txt.gz", "", basename(x))]
                                        return(z)
                                      }))


cyclo_2PG_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Cyclo2PG/", 
                                                 pattern = '*_norm.txt', full.names = T), 
                                      function(x) {
                                        z <- fread(x, col.names = 'avg_cvg')
                                        z[, sample := gsub("_norm.txt", "", basename(x))]
                                        return(z)
                                      }
))
cyclo_2PG_cvg <- merge(cyclo_2PG_raw_cvg, cyclo_2PG_cvg_nrm)








# -- coverage to Jm3101_v1.0
cyclo_HJmic_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jm3101_v1.0///", 
                                                  pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))

cyclo_HJmic_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jm3101_v1.0/", 
                                                      pattern = '*_norm.txt', full.names = T), 
                                           function(x) {
                                             z <- fread(x, col.names = 'avg_cvg')
                                             z[, sample := gsub("_norm.txt", "", basename(x))]
                                             return(z)
                                           }
))
cyclo_HJmic_cvg <- merge(cyclo_HJmic_raw_cvg, cyclo_HJmic_cvg_nrm)


# -- coverage to Jcali_primary
cyclo_HJcal_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jcali_primary/", 
                                                   pattern = '*.txt.gz', full.names = T),
                                        function(x){
                                          z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                          z[, sample := gsub(".txt.gz", "", basename(x))]
                                          return(z)
                                        }))

cyclo_HJcal_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jcali_primary//", 
                                                   pattern = '*_norm.txt', full.names = T), 
                                        function(x) {
                                          z <- fread(x, col.names = 'avg_cvg')
                                          z[, sample := gsub("_norm.txt", "", basename(x))]
                                          return(z)
                                        }
))
cyclo_HJcal_cvg <- merge(cyclo_HJcal_raw_cvg, cyclo_HJcal_cvg_nrm)







# -- coverage to Chandler
cyclo_hJ_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Chandler///", 
                                                pattern = '*.txt.gz', full.names = T),
                                     function(x){
                                       z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                       z[, sample := gsub(".txt.gz", "", basename(x))]
                                       return(z)
                                     }))

cyclo_hJ_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Chandler/", 
                                                pattern = '*_norm.txt', full.names = T), 
                                     function(x) {
                                       z <- fread(x, col.names = 'avg_cvg')
                                       z[, sample := gsub("_norm.txt", "", basename(x))]
                                       return(z)
                                     }
))
cyclo_hJ_cvg <- merge(cyclo_hJ_raw_cvg, cyclo_hJ_cvg_nrm)



# -- coverage to Pawnee
cyclo_HP_raw_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Pawnee_v1/", 
                                                pattern = '*.txt.gz', full.names = T),
                                     function(x){
                                       z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                       z[, sample := gsub(".txt.gz", "", basename(x))]
                                       return(z)
                                     }))

cyclo_HP_cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Pawnee_v1/", 
                                                pattern = '*_norm.txt', full.names = T), 
                                     function(x) {
                                       z <- fread(x, col.names = 'avg_cvg')
                                       z[, sample := gsub("_norm.txt", "", basename(x))]
                                       return(z)
                                     }
))
cyclo_HP_cvg <- merge(cyclo_HP_raw_cvg, cyclo_HP_cvg_nrm)





# ---- calculate coverage in 500bp windows

# ----- 2PA2023
C2PA_st <- 6000000
C2PA_en <- 7000000

cyclo_2PA_cvg[, window := cut(position, breaks = seq(C2PA_st, C2PA_en+1000, by = 1000), labels = seq(C2PA_st, C2PA_en, by = 1000), include.lowest =T), by = sample]
cyclo_2PA_cvg_1kb <- cyclo_2PA_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_2PA_cvg_1kb[, window := as.numeric(as.character((window)))]
cyclo_2PA_cvg_1kb[, nrm_cvg := coverage/avg_cvg]

cyclo_2PA_cvg_1kb <- merge(cyclo_2PA_cvg_1kb, samples, by = 'sample')


# ----- 2PG2023
C2PG_st <- 6500000
C2PG_en <- 6700000

cyclo_2PG_cvg[, window := cut(position, breaks = seq(C2PG_st, C2PG_en+1000, by = 1000), labels = seq(C2PG_st, C2PG_en, by = 1000), include.lowest =T), by = sample]
cyclo_2PG_cvg_1kb <- cyclo_2PG_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_2PG_cvg_1kb[, window := as.numeric(as.character((window)))]
cyclo_2PG_cvg_1kb[, nrm_cvg := coverage/avg_cvg]

cyclo_2PG_cvg_1kb <- merge(cyclo_2PG_cvg_1kb, samples, by = 'sample')





# ----- HJm3101
HJmic_st <- 30537000
HJmic_en <- 30588000

cyclo_HJmic_cvg[, window := cut(position, breaks = seq(HJmic_st, HJmic_en+500, by = 500), labels = seq(HJmic_st, HJmic_en, by = 500), include.lowest =T), by = sample]
cyclo_HJmic_cvg_500 <- cyclo_HJmic_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_HJmic_cvg_500[, window := as.numeric(as.character((window)))]
cyclo_HJmic_cvg_500[, nrm_cvg := coverage/avg_cvg]

cyclo_HJmic_cvg_500 <- merge(cyclo_HJmic_cvg_500, samples, by = 'sample')


# ----- HJcali_primary
HJcal_st <- 31348000
HJcal_en <- 31400000

cyclo_HJcal_cvg[, window := cut(position, breaks = seq(HJcal_st, HJcal_en+500, by = 500), labels = seq(HJcal_st, HJcal_en, by = 500), include.lowest =T), by = sample]
cyclo_HJcal_cvg_500 <- cyclo_HJcal_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_HJcal_cvg_500[, window := as.numeric(as.character((window)))]
cyclo_HJcal_cvg_500[, nrm_cvg := coverage/avg_cvg]

cyclo_HJcal_cvg_500 <- merge(cyclo_HJcal_cvg_500, samples, by = 'sample')



# ----- Chandler
cyclo_hJ_cvg <- cyclo_hJ_cvg[position >=31.86e6 & position <= 31.905e6]
hJ_st <- 31.86e6
hJ_en <- 31.905e6

cyclo_hJ_cvg[, window := cut(position, breaks = seq(hJ_st, hJ_en+500, by = 500), labels = seq(hJ_st, hJ_en, by = 500), include.lowest =T), by = sample]
cyclo_hJ_cvg_500 <- cyclo_hJ_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_hJ_cvg_500[, window := as.numeric(as.character((window)))]
cyclo_hJ_cvg_500[, nrm_cvg := coverage/avg_cvg]

cyclo_hJ_cvg_500 <- merge(cyclo_hJ_cvg_500, samples, by = 'sample')

# ----- HP
HP_st <- 6400000
HP_en <- 6800000

cyclo_HP_cvg[, window := cut(position, breaks = seq(HP_st, HP_en+500, by = 500), labels = seq(HP_st, HP_en, by = 500), include.lowest =T), by = sample]
cyclo_HP_cvg_500 <- cyclo_HP_cvg[, .(coverage = mean(coverage)), by = .(sample,window,avg_cvg)]
cyclo_HP_cvg_500[, window := as.numeric(as.character((window)))]
cyclo_HP_cvg_500[, nrm_cvg := coverage/avg_cvg]

cyclo_HP_cvg_500 <- merge(cyclo_HP_cvg_500, samples, by = 'sample')


# ===== plot

# ----- Cpal 2PA
ggplot(cyclo_2PA_cvg_1kb[ploidy == 'diploid'],
       aes(x = window, y = coverage, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
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
  )





# ----- Cpal 2PG
ggplot(cyclo_2PG_cvg_1kb[ ploidy == 'diploid'],
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(limits = c(0, 5)) +
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
  )








# ----- HJmic
#cyclo_HJ_cvg_500[sample == 'SRR23378890', phenotype := 'protogynous']
ggplot(cyclo_HJmic_cvg_500[ ploidy == 'diploid'],
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
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
  )


cyclo_HJ1mic_cvg_500 <- cyclo_HJmic_cvg_500[window > 30560000 & window < 30570000]
cyclo_HJ1mic_cvg_500[, nrm_cvg := coverage/avg_cvg]

ggplot(cyclo_HJ1mic_cvg_500,
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
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
  )

# ----- HJcali_primary
ggplot(cyclo_HJcal_cvg_500,
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
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
  )










# ----- hJ (Chandler)
ggplot(cyclo_hJ_cvg_500[ploidy == 'diploid'],
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized depth', color = '', title = '') +
  theme(aspect.ratio = 0.3,
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


# ----- HP
ggplot(cyclo_HP_cvg_500[window >6570000 & window < 6650000 & ploidy == 'diploid'],
       aes(x = window, y = nrm_cvg, group = sample, color = phenotype)) +
  geom_line() +
  #facet_grid(species~reference, scales = 'free_x') + 
  #geom_line(data = cvg500bp[genotype == '??'], linewidth = 0.5, alpha = 0.9)  +
  #geom_line(data = cvg500bp[genotype != '??'], linewidth = 0.5, alpha = 0.9)  +
  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  #scale_color_manual( values = c('tan', 'maroon', 'darkblue')) +
  
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
  )


