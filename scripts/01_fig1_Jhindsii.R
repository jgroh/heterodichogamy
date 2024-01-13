library(data.table)
library(ggplot2)
library(cowplot)
library(ggbreak)

# ----- read phenotypes -----
pheno <- fread("~/workspace/heterodichogamy/data/phenotypes.txt", header = F, col.names = c("ID", "phenotype", "species"))
pheno[phenotype == 'protogynous', genotype := 'H?']
pheno[phenotype == 'protandrous', genotype := 'hh']

# ----- read GWAS ------
x <- fread("~/workspace/heterodichogamy/hindsii/output/Putah2JcaliP.assoc.txt")

# snps_per_chr <- x[, .N, by = chr]
# setkey(snps_per_chr, N)
# chrs <- snps_per_chr[N > 10000, chr]
# 
chrs <- c('JAKSXK010000002.1', #1
          'JAKSXK010000003.1',
          'JAKSXK010000006.1',
          'JAKSXK010000009.1',
          'JAKSXK010000015.1',
          'JAKSXK010000005.1', #6
          'JAKSXK010000001.1',
          'JAKSXK010000010.1',
          'JAKSXK010000014.1',
          'JAKSXK010000008.1',
          'JAKSXK010000007.1',
          'JAKSXK010000011.1',
          'JAKSXK010000004.1', #13
          'JAKSXK010000013.1',
          'JAKSXK010000016.1',
          'JAKSXK010000012.1'
)
          
x <- x[chr %in% chrs]
x[, chr := factor(chr, levels=chrs)]
setkey(x, chr)
                                  
#x <-  x[!is.infinite(-log10(p_wald))]                               
x[, N:= seq_len(.N)]



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

hindsii_cvg_p <- merge(hindsii2cali_p_cvg, hindsii2cali_p_cvg_nrm)
hindsii_cvg_alt <- merge(hindsii2cali_alt_cvg, hindsii2cali_alt_cvg_nrm)
hindsii_cvg <- rbind(hindsii_cvg_p,hindsii_cvg_alt )








# ---- calculate coverage in 500bp windows -----
Jcal_pstart <- 31340000
Jcal_pend <- 31400000
Jcal_astart <- 30490000
Jcal_aend <- 30536000

winsize <- 1e3
hindsii_cvg[reference == 'primary', window := cut(position, breaks = seq(Jcal_pstart, Jcal_pend+winsize, by = winsize), labels = seq(Jcal_pstart, Jcal_pend, by = winsize), include.lowest =T), by = sample]
hindsii_cvg[reference == 'alternate', window := cut(position, breaks = seq(Jcal_astart, Jcal_aend+winsize, by = winsize), labels = seq(Jcal_astart, Jcal_aend, by = winsize), include.lowest =T), by = sample]

cvg_win <- hindsii_cvg[, .(coverage = mean(coverage)), by = .(sample, species, window, reference, avg_cvg)]
cvg_win[, window := as.numeric(as.character((window)))]

cvg_win <- merge(cvg_win, pheno[, .(sample = ID, phenotype, genotype)], all.x = T)
cvg_win[, nrm_cvg := coverage/avg_cvg]



# whole genome, Jcali primary assembly
wg_pvals <- ggplot(x, aes(x = N, y = -log10(p_lrt), color = chr)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = rep(c("black", "gray"), 8)) + 
  scale_y_continuous(limits = c(0, 210)) +
  scale_y_break(c(60, 190), ticklabels = c(190, 210)) +
  
  theme_classic() + 
  labs(y = expression(-log[10](P)), 
       x = '') +
  theme(
    #aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    axis.ticks.length.x = unit(0, 'cm'),
    axis.text.x = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(size = 10, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = -1),
    #legend.position = c(0.15, 0.85),
    legend.position = 'none'
  )
wg_pvals

ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/Jhindsii_gwas.pdf', plot = wg_pvals, width = 7, height = 1.8, units = 'in')


# H-locus, Jcali H assembly

pval_plt <- ggplot(x[chr == 'JAKSXK010000007.1' & ps > 31345000 & ps < 31400000], aes(x = ps, y = -log10(p_lrt))) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0, 210)) +
  scale_x_continuous(limits = c(31345000, 31400000), breaks = seq(31.35e6, 31.39e6, length.out = 3), labels = seq(31.35, 31.39, length.out = 3)) +
  theme_classic() +
  scale_y_break(c(60, 190), ticklabels = c(190, 210)) +
  labs(y = expression(-log[10](P)),
       x = '')  +
   annotate('rect', xmin = 31380614, xmax = 31383466,
            ymin = 55, ymax = 65, fill = 'blue', alpha = 0.2) +
   geom_segment(aes(x = 31383466, y = 59, xend = 31380614, yend = 59),
                arrow = arrow(length = unit(0.2, "cm")), color = 'blue') +
   annotate('rect', xmin = 31358581, xmax = 31360135,
            ymin = 55, ymax = 65, fill = 'maroon', alpha = 0.2) +
   geom_segment(aes(x = 31358581, y = 59, xend = 31360135, yend = 59),
                arrow = arrow(length = unit(0.2, "cm")), color = 'maroon') + 
  theme(
    #aspect.ratio = .3,
    plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(size = 12, vjust = -2),
    axis.title.y = element_text(size = 12, vjust = -1),
    #legend.position = c(0.15, 0.85),
    legend.position = 'none'
  )

pval_plt
ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/Jhindsii_pvals.pdf', plot = pval_plt, width = 7, height = 2, units = 'in')


# H-locus, Jcali alt assembly
# pval_plt <- ggplot(x[chr == 'JAKSXL010000006.1' & ps > 30.49e6 & ps < 30.529e6], aes(x = ps, y = -log10(p_wald))) + 
#   geom_point(size = 1) + 
#   scale_y_continuous(limits = c(0, 55)) +
#   theme_classic() + 
#   labs(y = expression(-log[10](P)), 
#        x = '') +
#   theme(aspect.ratio = .3,
#         plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12), 
#         axis.text = element_text(size = 12),
#         axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85), 
#         legend.position = 'none'
#   ) + 
#   annotate('rect', xmin = 30517523, xmax = 30520347, 
#            ymin = 40, ymax = 50, fill = 'blue', alpha = 0.2) +
#   geom_segment(aes(x = 30520347, y = 45, xend = 30517523, yend = 45),
#                arrow = arrow(length = unit(0.2, "cm")), color = 'blue') +
#   annotate('rect', xmin = 30498727, xmax = 30500052, 
#            ymin = 40, ymax = 50, fill = 'maroon', alpha = 0.2) + 
#   geom_segment(aes(x = 30498727, y = 45, xend = 30500052, yend = 45),
#                arrow = arrow(length = unit(0.2, "cm")), color = 'maroon') 


cvg_win[genotype == 'hh', genotype := 'gg']
cvg_win[genotype == 'H?', genotype := 'G?']
cvg_win[, genotype := factor(genotype, levels = c("gg", 'G?'))]


cvg_plt <- ggplot(cvg_win[genotype != '??' & species == 'hindsii' & reference == 'primary' & window > 313400000 & window < 31400000],
       aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  geom_line(data = cvg_win[ genotype != '??' & species == 'hindsii' & reference == 'primary' & window > 31340000 & window < 31400000], linewidth = 0.8, alpha = 0.9)  +

  #scale_color_manual( values = c('gray', 'maroon', 'tan', 'darkblue')) +
  scale_color_manual( values = c('tan', 'maroon')) +
  scale_x_continuous(limits = c(31345000, 31400000), breaks = seq(31.35e6, 31.39e6, length.out = 3), labels = seq(31.35, 31.39, length.out = 3)) +
  
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
    legend.key.size = unit(.5, 'cm')
  )

cvg_plt
ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/Jhindsii_cvg_H.pdf', plot = cvg_plt, width = 7, height = 2.5, units = 'in')




