library(data.table)
library(ggplot2)
library(cowplot)

# ----- read ids
sra <- fread("~/workspace/heterodichogamy/regia/sra_samples.txt")

# --- read coverage
regia_coverage <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Chandler/sra/", 
                                              pattern = "*.txt.gz",
                                              full.names = T),
       function(x){
         z <- fread(x, select = 2:3, col.names = c("pos", "coverage"))
         z[, sample := gsub(".txt.gz", "", basename(x))]
         return(z)
       }))
regia_coverage[, species := 'regia']


cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/regia/results/coverage_Chandler/sra/", 
                                       pattern = "*_norm.txt",
                                       full.names = T),
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }))
coverage <- merge(cvg_nrm, regia_coverage)


# ---- calculate coverage in 500bp windows
coverage[, window := cut(pos, breaks = seq(31.8e6, 32e6, by = 500), labels = seq(31.8e6, 32e6-500, by = 500), include.lowest =T), by = sample]
coverage[, window := as.numeric(as.character((window)))]
coverage_bins <- coverage[, .(coverage = mean(coverage)), by = .(sample, species, window, avg_cvg)]
coverage_bins[, nrm_cvg := coverage/avg_cvg]


# assign genotype based on CNV
#coverage_bins[window >= 31.865e6 & window <= 31.895e6 ]
sra_Hh <- coverage_bins[window >= 31.883e6 & window <= 31.885e6 & nrm_cvg > 2, unique(sample)]
sra_HH <- coverage_bins[window >= 31.883e6 & window <= 31.885e6 & nrm_cvg > 5, unique(sample)]
sra_Hh <- sra_Hh[!sra_Hh %in% sra_HH]
sra_hh <- coverage_bins[!sample %in% c(sra_Hh, sra_HH), unique(sample)]

coverage_bins[sample %in% sra_HH, genotype := 'HH']
coverage_bins[sample %in% sra_Hh, genotype := 'Hh']
coverage_bins[sample %in% sra_hh, genotype := 'hh']

# plot 15 at a time to verify it worked

coverage_bins[sample %in% coverage_bins[, unique(sample)][1:15], grp := 'grp1']
coverage_bins[sample %in% coverage_bins[, unique(sample)][16:30], grp := 'grp2']
coverage_bins[sample %in% coverage_bins[, unique(sample)][31:45], grp := 'grp3']
coverage_bins[sample %in% coverage_bins[, unique(sample)][46:60], grp := 'grp4']
coverage_bins[sample %in% coverage_bins[, unique(sample)][61:75], grp := 'grp5']
coverage_bins[sample %in% coverage_bins[, unique(sample)][76:87], grp := 'grp6']

ggplot(coverage_bins[window >= 31.865e6 & window <= 31.895e6 ],
       aes(x = window, y = nrm_cvg, group = sample, color = genotype)) +
  geom_line(linewidth = 0.8) +
  #geom_hline(yintercept = 5) +
  
  theme_classic() + 
  labs(x = 'Position (Mb)', y =  'Normalized read depth') +
  theme(aspect.ratio = 0.25,
        legend.position = c(0.1, 0.9),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size = 14, vjust = -2),
        legend.text = element_text(size = 12)
  )



# ----- merge sra sample list with genotypes

sra_geno <- merge(sra, coverage_bins[window == min(window), .(run = sample, genotype)])

fwrite(sra_geno, file = '~/workspace/heterodichogamy/regia/sra_samples_geno.txt', row.names = F, quote = F, sep = "\t", col.names = T)


