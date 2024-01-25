
library(ggplot2)
library(cowplot)
library(data.table)


# ----  Cyclocarya coverage to Jcali primary -----
Cpal_meta <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/wgs_samples.tsv")

Cpal_JugG <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jcali_primary/", 
                                           pattern = '*.txt.gz', full.names = T),
                                function(x){
                                  z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                  z[, sample := gsub(".txt.gz", "", basename(x))]
                                  return(z)
                                }))

Cpal_JugG_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Cyclocarya_paliurus/results/coverage_Jcali_primary/", 
                                           pattern = '*_norm.txt', full.names = T), 
                                function(x) {
                                  z <- fread(x, col.names = 'avg_cvg')
                                  z[, sample := gsub("_norm.txt", "", basename(x))]
                                  return(z)
                                }
))

Cpal_cvg <- merge(Cpal_JugG, Cpal_JugG_nrm)
Cpal_cvg <- Cpal_cvg[sample %in% Cpal_meta[ploidy == 'diploid', run]]
Cpal_cvg <- merge(Cpal_cvg, Cpal_meta[, .(sample = run, ID)])

Cpal_cvg[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Cyc <- Cpal_cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg, ID)]
cvg1kb_Cyc[, window := as.numeric(as.character((window)))]
cvg1kb_Cyc[, nrm_cvg := coverage/avg_cvg]


Cpal_plot <- ggplot(cvg1kb_Cyc,
       aes(x = window, y = nrm_cvg, group = ID, color = ID)) +
  geom_line() +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized\nread depth', color = '', title = 'Cyclocarya paliurus') +
  theme(aspect.ratio = 0.25,
        plot.title = element_text(face = 'italic'),
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
        axis.title.y = element_text(size=12),
        legend.position = 'none'

  ) + 
  annotate("rect", xmin = 31380787, xmax = 31383637, ymin = 3.5, ymax = 4,
         alpha = .2,fill = 'blue') + 
  geom_segment(aes(x = 31383637, y = 3.75, xend = 31380787, yend = 3.75),
               arrow = arrow(length = unit(0.15, "cm")), color = 'blue', linewidth = 0.2) 




# ------- Pterocarya ----- 


Pste_JugG <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Pterocarya_stenoptera/results/coverage_Jcali_primary/", 
                                         pattern = '*.txt.gz', full.names = T),
                              function(x){
                                z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                z[, sample := gsub(".txt.gz", "", basename(x))]
                                return(z)
                              }))

Pste_JugG_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Pterocarya_stenoptera/results/coverage_Jcali_primary/", 
                                             pattern = '*_norm.txt', full.names = T), 
                                  function(x) {
                                    z <- fread(x, col.names = 'avg_cvg')
                                    z[, sample := gsub("_norm.txt", "", basename(x))]
                                    return(z)
                                  }
))

Pste_cvg <- merge(Pste_JugG, Pste_JugG_nrm)

Pste_cvg[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Pste <- Pste_cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Pste[, window := as.numeric(as.character((window)))]
cvg1kb_Pste[, nrm_cvg := coverage/avg_cvg]


Pste_plot <- ggplot(cvg1kb_Pste,
       aes(x = window, y = nrm_cvg, group = sample, color = sample)) +
  geom_line() +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized\nread depth', color = '', title = 'Pterocarya stenoptera') +
  theme(aspect.ratio = 0.25,

        plot.title = element_text(face = 'italic'),
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
        axis.title.y = element_text(size=12),
        legend.position = 'none'
  ) + 
  annotate("rect", xmin = 31380787, xmax = 31383637, ymin = 3.5, ymax = 4,
           alpha = .2,fill = 'blue') + 
  geom_segment(aes(x = 31383637, y = 3.75, xend = 31380787, yend = 3.75),
               arrow = arrow(length = unit(0.15, "cm")), color = 'blue', linewidth = 0.2) 



# ------- Platycarya ----- 


Pstr_JugG <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Platycarya_strobilaceae/results/coverage_Jcali_primary/", 
                                         pattern = '*.txt.gz', full.names = T),
                              function(x){
                                z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                z[, sample := gsub(".txt.gz", "", basename(x))]
                                return(z)
                              }))

Pstr_JugG_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Platycarya_strobilaceae/results/coverage_Jcali_primary/", 
                                             pattern = '*_norm.txt', full.names = T), 
                                  function(x) {
                                    z <- fread(x, col.names = 'avg_cvg')
                                    z[, sample := gsub("_norm.txt", "", basename(x))]
                                    return(z)
                                  }
))

Pstr_cvg <- merge(Pstr_JugG, Pstr_JugG_nrm)

Pstr_cvg[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Pstr <- Pstr_cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Pstr[, window := as.numeric(as.character((window)))]
cvg1kb_Pstr[, nrm_cvg := coverage/avg_cvg]


Pstr_plot <- ggplot(cvg1kb_Pstr,
                    aes(x = window, y = nrm_cvg, group = sample, color = sample)) +
  geom_line() +
  theme_classic() + 
  scale_x_continuous(n.breaks = 5) +
  labs(x = '', y =  'Normalized\nread depth', color = '', title = 'Platycarya strobilaceae') +
  theme(aspect.ratio = 0.25,
        
        plot.title = element_text(face = 'italic'),
        axis.text.x = element_text(size=10, angle = 60, vjust = 0.6),
        axis.title.y = element_text(size=12),
        legend.position = 'none'
  ) + 
  annotate("rect", xmin = 31380787, xmax = 31383637, ymin = 3.5, ymax = 4,
           alpha = .2,fill = 'blue') + 
  geom_segment(aes(x = 31383637, y = 3.75, xend = 31380787, yend = 3.75),
               arrow = arrow(length = unit(0.15, "cm")), color = 'blue', linewidth = 0.2) 








plot_grid(Cpal_plot, Pste_plot, Pstr_plot, ncol = 1, align = 'v')





