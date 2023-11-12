library(data.table)
library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)


# ---- read GWAS data -----
gbs <- fread("~/workspace/heterodichogamy/pecan/GBS_cultivars_gemma.assoc.txt")
gbs[, N:= seq_len(.N)]


# -----  read phenotypes of WGS data
fam <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.fam")
fam <- fam[, .(V1, V6)]
setnames(fam, c("V1", "V6"), c("sample", "phenotype"))

wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_reassigned.tsv")
setnames(wgs_samples, 'run', 'sample')

# ----- read LD data -----
ld <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.geno.ld")
setnames(ld, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))

# ----- read HC locus genes
gn <- fread("~/workspace/heterodichogamy/HC_locus_structure/HC_locus_genes.txt")

# -------- library IDs of 18 sp

sra <- fread("~/workspace/heterodichogamy/Carya_15sp/Carya_15sp_libraries.txt")

# --  read coverage of WGS data
wgs_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Pawnee/", 
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Pawnee//", 
                                       pattern = '*_norm.txt', full.names = T), 
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }
))

cvg <- merge(wgs_cvg, cvg_nrm)
cvg <- merge(cvg, wgs_samples)

# ------ read coverage of 15 sp ------
cvg15sp <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Carya_15sp/results/coverage_Lakota_v1//", 
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))
cvg15sp_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Carya_15sp/results/coverage_Lakota_v1//", 
                                       pattern = '*_norm.txt', full.names = T), 
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }
))
cvg15sp <- merge(cvg15sp, cvg15sp_nrm)


# ----- Plot GWAS

# whole genome 
gwas_wg <- ggplot(gbs, aes(x = N, y = -log10(p_wald), color = chr)) + 
         geom_point(size =1) + 
         scale_color_manual(values = rep(c("tan", "purple4"), 8)) + 
  labs(y = expression(-log[10](P)), 
       x = '') +
  theme_classic() + 
  theme(aspect.ratio = .3,
       # plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2),
        #legend.position = c(0.15, 0.85), 
        legend.position = 'none')



# ----- zoomed in on GWAS peak -----
gwas_peak <- ggplot(gbs[chr == 'CM031812.1' & ps > 6.4e6 & ps < 6.75e6], aes(x = ps, y = -log10(p_wald))) + 
  geom_point() + 
  #scale_color_manual(values = rep(c("black", "grey"), 8)) + 
  theme_classic() +
  labs(y = expression(-log[10](P)), 
       x = '') +
  theme(aspect.ratio = .3,
      #plot.margin = margin(30,30,30,30, "pt"),
      text = element_text(size = 12), 
      axis.text = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 12, vjust = -2),
      axis.title.y = element_text(size = 12, vjust = 2),
      #legend.position = c(0.15, 0.85), 
      legend.position = 'none')
  #annotate("rect", xmin = 6534137, xmax = 6534830, ymin = 0, ymax = 25,
  #          alpha = .8,fill = '#94435b') + 
  #annotate("rect", xmin = 6632000, xmax = 6634000, ymin = 0, ymax = 25,
  #         alpha = .8,fill = '#94435b') 




# ---- calculate coverage in 1000bp windows
st <- 6000000
en <- 7000000
cvg[, window := cut(position, breaks = seq(st, en+1000, by = 1000), labels = seq(st, en, by = 1000), include.lowest =T), by = sample]

cvg1kb <- cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg, type, phenotype)]
cvg1kb[, window := as.numeric(as.character((window)))]

cvg1kb[, nrm_cvg := coverage/avg_cvg]

# looks like there presence/absence variation
ggplot(cvg1kb[window >6400000 & window < 6700000], 
       aes(x = window, y = nrm_cvg, group = sample)) + geom_line()



# ---- -calculate coverage for 15sp
cvg15sp[, window := cut(position, breaks = seq(st, en+1000, by = 1000), labels = seq(st, en, by = 1000), include.lowest =T), by = sample]

cvg15sp1kb <- cvg15sp[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg15sp1kb[, window := as.numeric(as.character((window)))]

cvg15sp1kb[, nrm_cvg := coverage/avg_cvg]
cvg15sp1kb <- merge(cvg15sp1kb, sra[, .(sample = Run, Species)])

ggplot(cvg15sp1kb[window > 6500000 & window <7000000 & Species == 'ovata'], 
       aes(x = window, y = nrm_cvg, 
           group = sample)) + geom_line() + scale_y_continuous(limits = c(0,2))

# ----- color by phenotype *** note that there is clear evidence that some phenotypes were mismatched in the wgs data set ****
#ggplot(cvg1kb[window >6640000 & window < 6642500], 
#       aes(x = window, y = nrm_cvg, 
#           group = sample, color = as.factor(phenotype))) + geom_line()

# which samples don't match???
# cvg1kb[window >6640000 & window < 6642500 & nrm_cvg < 0.8 & nrm_cvg > 0.2, unique(type)]

# ----- color by phenotype *** FOR SUBSET OF FILTERED INDIVIDUALS BASED ON RELATEDNESS TO GBS SAMPLES ****

# type %in% c("Barton", "Dependable", "Excel", "Farley","Kiowa", "Mahan", "Oconee", "Pawnee", "Shoshoni","Sioux", "Sumner", "VC1-68"
ggplot(cvg1kb[window >6640000 & window < 6642500], 
       aes(x = window, y = nrm_cvg, 
           group = sample, color = as.factor(phenotype))) + geom_line()

#highconf_samples <- c("Barton", "Dependable", "Excel", "Farley",
#                      "Kiowa", "Mahan", "Oconee", "Pawnee", "Shoshoni",
#                      "Sioux", "Sumner", "VC1-68")
#cvg1kb[!type %in% highconf_samples, genotype := "NA"]
cvg1kb[type == 'Mahan', genotype := 'HH']
cvg1kb[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
cvg1kb[phenotype == 'protandrous', genotype := 'hh']

# ----- assign SV SV genotype *** for all samples, including those which phenotypes are off for 
ggplot(cvg1kb[window >6640000 & window < 6642500], 
       aes(x = window, y = nrm_cvg, 
           group = sample)) + geom_line()
SV_HH <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg < 0.1, unique(sample)]
cvg1kb[sample %in% SV_HH, SV := 2]

SV_Hh <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg > 0.1 & nrm_cvg < 0.8, unique(sample)]
cvg1kb[sample %in% SV_Hh, SV := 1]

SV_hh <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg > 0.8, unique(sample)]
cvg1kb[sample %in% SV_hh, SV := 0]

# check that the SV is scored correctly
ggplot(cvg1kb[window >6640000 & window < 6642500], 
       aes(x = window, y = nrm_cvg, 
           group = sample, color = as.factor(SV))) + geom_line()


# now zoom out again
# plot3 <- ggplot(cvg1kb[window >6400000 & window < 6750000], 
#             aes(x = window, y = nrm_cvg, group = sample, 
#                 color = as.factor(SV))) + geom_line() + 
#   scale_y_continuous(limits = c(0, 3)) +
#   scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
#   theme_classic() + 
#   labs(y = 'Normalized read depth', x = '') +
#   theme(
#     aspect.ratio = .3,
#     #plot.margin = margin(30,30,30,30, "pt"),
#     text = element_text(size = 12), 
#     axis.text = element_text(size = 12),
#     #axis.text.x = element_blank(),
#     axis.title.x = element_text(size = 12, vjust = -2),
#     axis.title.y = element_text(size = 12, vjust = 2),
#     #legend.position = c(0.15, 0.85), 
#     legend.position = 'none')

# ld[pos1 > 6.45e6 & pos1 < 6.67e6 & pos2 >  6.45e6 & pos2 < 6.67e6], # full region

ld_plt <- ggplot(ld[pos1 > 6.45e6 & pos1 < 6.67e6 & pos2 >  6.45e6 & pos2 < 6.67e6], 
                 aes(x = pos2, y = pos1, fill = r2)) + 
  geom_tile(width = 200, height = 200) +
  #annotate("rect", xmin = 31868751, xmax = 31870103, 
  #         ymin = 31868751, ymax = 31870103,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  #annotate("rect", xmin = 31884268, xmax = 31887072, 
  #         ymin = 31884268, ymax = 31887072,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "", y = "Pawnee chr 4 (Mb)", fill = expression(r^2)) +
  scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 9), labels = seq(6.45, 6.65, length.out = 9)) +
  scale_y_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 9), labels = seq(6.45, 6.65, length.out = 9)) +
  theme_classic() + 
  geom_hline(yintercept = 6462000, color ='red') +
  geom_hline(yintercept = 6662000, color = 'red') +
  geom_vline(xintercept = 6462000, color ='red') +
  geom_vline(xintercept = 6662000, color = 'red') +
  theme(aspect.ratio = 0.3, 
        legend.position = c(0.1, 0.6), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        plot.margin = unit(c(0,0,0,2), 'cm'))

ld_plt

cvg_plt <- ggplot(cvg1kb[ window >6450000 & window < 6670000],
            ) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample,color = genotype)) +
  #scale_y_continuous(limits = c(0, 3)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = seq(6.45, 6.65, length.out = 5))  +
  geom_vline(aes(xintercept = 6.662e6)) +
  theme_classic() +
  labs(y = 'Normalized read depth', x = 'Pawnee chr 4 (Mb)') +
  theme(
    aspect.ratio = .3,
    #plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, vjust = 2),
    #legend.position = c(0.15, 0.85),
    legend.position = 'none', 
    plot.margin = unit(c(0,0,0,2), 'cm')) + 
  geom_rect(data = gn, aes(xmin = Start, xmax = Stop, ymin=1, ymax=1.5), alpha = 0.3) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample,color = genotype)) 
  


#cvg_plt

#grid.arrange(ld_plt, cvg_plt, ncol = 1)
plot_grid(ld_plt, cvg_plt, ncol=1, align='v')

ggarrange(ld_plt, cvg_plt, ncol = 1)



cvg_plt2 <- ggplot(cvg15sp1kb[window >6450000 & window < 6670000]) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample), linewidth = 0.1) +
  scale_y_continuous(limits = c(0, 3)) +
  #scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = seq(6.45, 6.65, length.out = 5))  +
  theme_classic() +
  labs(y = 'Normalized read depth', x = 'Pawnee chr 4 (Mb)') +
  theme(
    aspect.ratio = .3,
    #plot.margin = margin(30,30,30,30, "pt"),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    #axis.text.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, vjust = 2),
    #legend.position = c(0.15, 0.85),
    legend.position = 'none', 
    plot.margin = unit(c(0,0,0,2), 'cm')) + 
  geom_rect(data = gn, aes(xmin = Start, xmax = Stop, ymin=1, ymax=1.5), alpha = 0.3) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample)) 

cvg_plt2
cvg_plt2

plot_grid(cvg_plt, cvg_plt2, ncol = 1, align='v')


# test for association between SV and phenotype
testdata <- cvg1kb[type %in% highconf_samples, .SD[1], by = type]
testdata[, SV_present := ifelse(SV > 0, 1, 0)]
chisq.test(x = testdata$SV_present, y = testdata$phenotype)

plot1
plot2
plot3



ggarrange(plot1, plot2, plot3, ncol = 1)

