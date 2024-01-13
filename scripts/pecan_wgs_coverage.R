
library(ggplot2)
library(cowplot)
library(data.table)


# ----- 
wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c('sample', 'genotype'))


# ----  coverage to Pawnee -----
wgs_cvg_Paw <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Pawnee/", 
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

cvg_nrm_Paw <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Pawnee/", 
                                       pattern = '*_norm.txt', full.names = T), 
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }
))

cvg_Paw <- merge(wgs_cvg_Paw, cvg_nrm_Paw)




# ------ coverage to Lakota v1
wgs_cvg_Lak1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_v1/", 
                                                  pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))
cvg_nrm_Lak1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_v1/", 
                                                      pattern = '*_norm.txt', full.names = T), 
                                           function(x) {
                                             z <- fread(x, col.names = 'avg_cvg')
                                             z[, sample := gsub("_norm.txt", "", basename(x))]
                                             return(z)
                                           }
))
cvg_Lak1 <- merge(wgs_cvg_Lak1, cvg_nrm_Lak1)

# ------ coverage to Lakota alt
wgs_cvg_Lak2 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_alt/", 
                                            pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
cvg_nrm_Lak2 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_alt/", 
                                            pattern = '*_norm.txt', full.names = T), 
                                 function(x) {
                                   z <- fread(x, col.names = 'avg_cvg')
                                   z[, sample := gsub("_norm.txt", "", basename(x))]
                                   return(z)
                                 }
))
cvg_Lak2 <- merge(wgs_cvg_Lak2, cvg_nrm_Lak2)


# ------ coverage to Elliott v1
wgs_cvg_Ell1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Elliott_v1", 
                                            pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
cvg_nrm_Ell1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Elliott_v1/", 
                                            pattern = '*_norm.txt', full.names = T), 
                                 function(x) {
                                   z <- fread(x, col.names = 'avg_cvg')
                                   z[, sample := gsub("_norm.txt", "", basename(x))]
                                   return(z)
                                 }
))
cvg_Ell1 <- merge(wgs_cvg_Ell1, cvg_nrm_Ell1)


# ------ coverage to Oaxaca
wgs_cvg_Oax1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Oaxaca_v1/", 
                                            pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
cvg_nrm_Oax1 <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Oaxaca_v1//", 
                                            pattern = '*_norm.txt', full.names = T), 
                                 function(x) {
                                   z <- fread(x, col.names = 'avg_cvg')
                                   z[, sample := gsub("_norm.txt", "", basename(x))]
                                   return(z)
                                 }
))
cvg_Oax1 <- merge(wgs_cvg_Oax1, cvg_nrm_Oax1)


# ------ coverage to Csin
wgs_cvg_Csin <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Csin_BNU/", 
                                            pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
cvg_nrm_Csin <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Csin_BNU//", 
                                            pattern = '*_norm.txt', full.names = T), 
                                 function(x) {
                                   z <- fread(x, col.names = 'avg_cvg')
                                   z[, sample := gsub("_norm.txt", "", basename(x))]
                                   return(z)
                                 }
))
cvg_Csin <- merge(wgs_cvg_Csin, cvg_nrm_Csin)

# ------ coverage to Ccat
wgs_cvg_Ccat <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Ccat_BNU/", 
                                            pattern = '*.txt.gz', full.names = T),
                                 function(x){
                                   z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                   z[, sample := gsub(".txt.gz", "", basename(x))]
                                   return(z)
                                 }))
cvg_nrm_Ccat <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Ccat_BNU//", 
                                            pattern = '*_norm.txt', full.names = T), 
                                 function(x) {
                                   z <- fread(x, col.names = 'avg_cvg')
                                   z[, sample := gsub("_norm.txt", "", basename(x))]
                                   return(z)
                                 }
))
cvg_Ccat <- merge(wgs_cvg_Ccat, cvg_nrm_Ccat)




# ----- Pawnee -----
#6640000 & window < 6642500
cvg_Paw[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Paw <- cvg_Paw[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Paw[, window := as.numeric(as.character((window)))]
cvg1kb_Paw[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Paw <- merge(cvg1kb_Paw, wgs_samples)

cvg1kb_Paw[genotype == 'HH', genotype := 'GG']
cvg1kb_Paw[genotype == 'Hh', genotype := 'Gg']
cvg1kb_Paw[genotype == 'hh', genotype := 'gg']

pawnee_cvg <- ggplot(cvg1kb_Paw[window > 6.42e6 & window < 6.7e6], aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  geom_line() + 
  #geom_line(data = cvg1kb_Lak1[sample == 'SRR15911540' & window > 6400000 & window < 7100000], color = 'black' ) +
  scale_color_manual(values = c( 'tan', 'maroon', 'turquoise4')) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 3)) +
  labs(x = 'Chromosome 4', y = 'Normalized\nread depth', color = '', title = "'Pawnee'") + 
  geom_vline(xintercept = 6461560, linetype = 2) + 
  geom_vline(xintercept = 6658636, linetype = 2)  + 
  theme(aspect.ratio = .3,
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.9, 0.8),
        plot.title = element_text(size = 12),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(face = 'italic'))

pawnee_cvg






# ---- Lakota v1 -----

cvg_Lak1[, window := cut(position, breaks = seq(min(position), max(position)+5000, by = 5000), labels = seq(min(position), max(position), by = 5000), include.lowest =T), by = sample]
cvg1kb_Lak1 <- cvg_Lak1[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Lak1[, window := as.numeric(as.character((window)))]
cvg1kb_Lak1[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Lak1 <- merge(cvg1kb_Lak1, wgs_samples)
#cvg1kb_Lak1[type == 'Mahan', genotype := 'HH']
#cvg1kb_Lak1[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
#cvg1kb_Lak1[phenotype == 'protandrous', genotype := 'hh']

ggplot(cvg1kb_Lak1[window > 6400000 & window < 7100000], aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  geom_line() + 
  #geom_line(data = cvg1kb_Lak1[sample == 'SRR15911540' & window > 6400000 & window < 7100000], color = 'black' ) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  theme_classic() + 
  theme(aspect.ratio = .3) + 
  labs(x = 'Position', y = 'Normalized read depth', color = '') + 
  geom_vline(xintercept = 6507677, linetype = 2) + 
  geom_vline(xintercept = 6952981, linetype = 2) 
  
  


# ---- Lakota alt

cvg_Lak2[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Lak2 <- cvg_Lak2[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Lak2[, window := as.numeric(as.character((window)))]
cvg1kb_Lak2[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Lak2 <- merge(cvg1kb_Lak2, wgs_samples)

cvg1kb_Lak2[genotype == 'HH', genotype := 'GG']
cvg1kb_Lak2[genotype == 'Hh', genotype := 'Gg']
cvg1kb_Lak2[genotype == 'hh', genotype := 'gg']

lak2_cvg <- ggplot(cvg1kb_Lak2, aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  scale_color_manual(values = c('tan', 'maroon', 'turquoise4')) +
  scale_y_continuous(limits = c(0, 2)) +
  geom_line() + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  labs(x = 'Scaffold 104', y = 'Normalized\nread depth', title = "'Lakota' alternate assembly") + 
  theme(aspect.ratio = .3,
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, vjust = 2),
        legend.position = 'none',
        plot.title = element_text(size = 12),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(face = 'italic')) + 
  geom_vline(xintercept = 286632, linetype = 2)  +
  geom_vline(xintercept = 64838, linetype = 2) 


#286632-64838


plot_grid(pawnee_cvg, lak2_cvg, ncol = 1, align = 'v')






# ---- Elliott v1

cvg_Ell1[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Ell1 <- cvg_Ell1[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Ell1[, window := as.numeric(as.character((window)))]
cvg1kb_Ell1[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Ell1 <- merge(cvg1kb_Ell1, wgs_samples)
cvg1kb_Ell1[type == 'Mahan', genotype := 'HH']
cvg1kb_Ell1[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
cvg1kb_Ell1[phenotype == 'protandrous', genotype := 'hh']

ggplot(cvg1kb_Ell1, aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  geom_line()

ggplot(cvg1kb_Ell1[window > 6200000], aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  geom_line()


# ---- Oaxaca

cvg_Oax1[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Oax1 <- cvg_Oax1[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Oax1[, window := as.numeric(as.character((window)))]
cvg1kb_Oax1[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Oax1 <- merge(cvg1kb_Oax1, wgs_samples)


ggplot(cvg1kb_Oax1, aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  geom_line()

ggplot(cvg1kb_Ell1[window > 6200000], aes(x = window, 
                                          y = nrm_cvg, 
                                          group = sample, 
                                          color = genotype)) + 
  geom_line()


# -----------

cvg_Csin[, window := cut(position, breaks = seq(min(position), max(position)+5000, by = 5000), labels = seq(min(position), max(position), by = 5000), include.lowest =T), by = sample]
cvgXkb_Csin <- cvg_Csin[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvgXkb_Csin[, window := as.numeric(as.character((window)))]
cvgXkb_Csin[, nrm_cvg := coverage/avg_cvg]

cvgXkb_Csin <- merge(cvgXkb_Csin, wgs_samples)

cvgXkb_Csin[genotype == 'HH', genotype := 'GG']
cvgXkb_Csin[genotype == 'Hh', genotype := 'Gg']
cvgXkb_Csin[genotype == 'hh', genotype := 'gg']

Csin_plot <- ggplot(cvgXkb_Csin[window > 5.9e6 & window < 6.6e6], aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks = seq(6e6, 6.6e6, length.out = 4), labels = seq(6, 6.6, length.out = 4))  +
  scale_color_manual(values = c('tan','maroon', 'turquoise4')) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic() + 
  labs(x = 'Csin G hap chr4 (Mb)', y = 'Normalized\nread depth', color = '') +
  theme(aspect.ratio = 0.2, 
        legend.position = c(0.8, 1.1),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        axis.title.y = element_text(),
        axis.title.x = element_text(), 
        legend.text = element_text(face = 'italic'),
        legend.key.size = unit(0.4, 'cm')
  )  + 
  geom_vline(xintercept = 6074902, linetype = 2) + 
  geom_vline(xintercept = 6421582, linetype = 2)

Csin_plot

# ----------- cathayensis

cvg_Ccat[, window := cut(position, breaks = seq(min(position), max(position)+5000, by = 5000), labels = seq(min(position), max(position), by = 5000), include.lowest =T), by = sample]
cvgXkb_Ccat <- cvg_Ccat[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvgXkb_Ccat[, window := as.numeric(as.character((window)))]
cvgXkb_Ccat[, nrm_cvg := coverage/avg_cvg]

cvgXkb_Ccat <- merge(cvgXkb_Ccat, wgs_samples)
cvgXkb_Ccat[genotype == 'HH', genotype := 'GG']
cvgXkb_Ccat[genotype == 'Hh', genotype := 'Gg']
cvgXkb_Ccat[genotype == 'hh', genotype := 'gg']


Ccat_plot <- ggplot(cvgXkb_Ccat[window > 7.05e6 & window < 7.55e6], aes(x = window, 
                                                         y = nrm_cvg, 
                                                         group = sample, 
                                                         color = genotype)) + 
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks = seq(7.1e6, 7.5e6, length.out = 5), labels = seq(7.1, 7.5, length.out = 5))  +
  scale_color_manual(values = c('tan','maroon', 'turquoise4')) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic() + 
  labs(x = 'Ccat g hap chr4 (Mb)', y = 'Normalized\nread depth', color = '') +
  theme(aspect.ratio = 0.2, 
        legend.position = 'none',
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        axis.title.y = element_text(),
        axis.title.x = element_text()
  )  + 
  geom_vline(xintercept = 7186465, linetype = 2) + 
  geom_vline(xintercept = 7419255, linetype = 2)

Ccat_plot

#library(gridExtra)
#library(ggpubr)
plot_grid(Csin_plot, Ccat_plot, ncol = 1, align = 'v')
