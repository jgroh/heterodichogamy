library(data.table)
library(data.table)
library(ggplot2)
library(cowplot)


# ----- 
wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples.tsv")
setnames(wgs_samples, 'run', 'sample')


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



# ---- Lakota v1

cvg_Lak1[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Lak1 <- cvg_Lak1[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Lak1[, window := as.numeric(as.character((window)))]
cvg1kb_Lak1[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Lak1 <- merge(cvg1kb_Lak1, wgs_samples)
cvg1kb_Lak1[type == 'Mahan', genotype := 'HH']
cvg1kb_Lak1[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
cvg1kb_Lak1[phenotype == 'protandrous', genotype := 'hh']

ggplot(cvg1kb_Lak1, aes(x = window, y = nrm_cvg, group = sample, color = genotype)) + geom_line()


# ---- Lakota alt

cvg_Lak2[, window := cut(position, breaks = seq(min(position), max(position)+1000, by = 1000), labels = seq(min(position), max(position), by = 1000), include.lowest =T), by = sample]
cvg1kb_Lak2 <- cvg_Lak2[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
cvg1kb_Lak2[, window := as.numeric(as.character((window)))]
cvg1kb_Lak2[, nrm_cvg := coverage/avg_cvg]

cvg1kb_Lak2 <- merge(cvg1kb_Lak2, wgs_samples)
cvg1kb_Lak2[type == 'Mahan', genotype := 'HH']
cvg1kb_Lak2[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
cvg1kb_Lak2[phenotype == 'protandrous', genotype := 'hh']

ggplot(cvg1kb_Lak2, aes(x = window, 
                        y = nrm_cvg, 
                        group = sample, 
                        color = genotype)) + 
  #scale_y_continuous(limits = c(0, 2)) +
  geom_line()


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

