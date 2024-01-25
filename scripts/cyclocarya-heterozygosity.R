library(data.table)
library(ggplot2)
library(magrittr)



ids <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/calls/diploids_2PG_CM056950.1_filtered.012.indv", header = F)
pos <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/calls/diploids_2PG_CM056950.1_filtered.012.pos", header = F)
het0 <- fread("~/workspace/heterodichogamy/Cyclocarya_paliurus/calls/diploids_2PG_CM056950.1_filtered.012", header = F)


# reformat
het0[, V1 := ids$V1]
setnames(het0, c("ID", pos$V2))
het1 <- melt(het0, id.vars = 'ID', variable.name = "pos", value.name = 'genotype')
het1[, pos := as.numeric(as.character(pos))]
head(het1)

het1[, window := cut(pos, breaks = seq(min(pos), max(pos)+5000, by = 5000), labels = seq(min(pos), max(pos), by = 5000), include.lowest =T), by = ID]
het1[, het := ifelse(genotype == 1, 1, 0)]
het1[genotype == -1, het := NA]
het2 <- het1[, .(het1kb = mean(het, na.rm = T)), by = .(ID, window)]
het2[, window := as.numeric(as.character((window)))]


ggplot(het2[window > 6.5e6 & window < 6.9e6], aes(x = window, y = het1kb, color = ID)) + 
  geom_point() + 
  theme_classic() +  
 geom_smooth(aes(group = ID), method = 'loess', span = 0.2, se = F)

possible_hets <- c('SRR6804841', 'SRR6804842','SRR6804843', 'SRR6804844','SRR6804845', 'SRR6804847')
ggplot(het2.1[window > 6.3e6 & window < 6.8e6 & ID %in% possible_hets], aes(x = window, y = het1kb, color = ID)) + 
  geom_point() + 
  theme_classic() + 
  geom_smooth(aes(group = ID), method = 'loess', span = 0.4, se = F)


