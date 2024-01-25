# ====== Approach 1 ======
# estimate pi using just homozygotes. problem is low sample size for G haplotype

p <- fread("~/workspace/heterodichogamy/pecan/WGS_data/GG-vs_gg_Lakota1_pi.txt")
p[, .(pi = weighted.mean(x = avg_pi, w = window_pos_2 - window_pos_1, na.rm=T)), by = pop]


gg_pi_vec <- sapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/pixy_gg_resampling/", full.names = T), 
       function(x){
         z <- fread(x)
         return(z[, weighted.mean(x = avg_pi, w = window_pos_2 - window_pos_1, na.rm=T)])
       }
)
names(gg_pi_vec) <- NULL
any(gg_pi_vec < 0.001812456)

1/78

#coords <- fread("~/workspace/heterodichogamy/pecan/WGS_data/Lakota_v1_Gloc_CDS_sites_for_haplotype_pi.txt")
#coords

# ===== Approach 2 =====
# estimate pi from phased haplotypes 
#p2 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/pixy_Gloc_CDS_haplotypes_pi.txt")

#p2[, mean(avg_pi), by = pop]




#p <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/WGS_cultivars_HCloc_CDS_pi.txt")
#d <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/WGS_cultivars_homozygotes_HCloc_CDS_dxy.txt")
#d2 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/WGS_cultivars_Hh_vs_hh_HCloc_CDS_dxy.txt")

#plot_data <- rbind(d[, .(pi = mean(avg_dxy), pop = 'HH_vs_hh'), by = .(window_pos_1, window_pos_2, count_comparisons)],
#      p[, .(pi = mean(avg_pi)), by = .(window_pos_1, window_pos_2, pop, count_comparisons)])
      
#ggplot(plot_data, aes(x = window_pos_2 - window_pos_1, y = count_comparisons)) + geom_point()
#ggplot(plot_data, aes(x = pop, y = pi, group = pop, color = pop)) + 
#  geom_point() + 
#  geom_boxplot()

#ggplot(plot_data[pop %in% c("hh", 'HH')], aes(x = pop, y = pi, group = pop, color = pop)) + 
#  geom_point() + 
#  geom_boxplot()

#plot_data[pop %in% c("hh", 'HH') & !is.na(pi), .N, by = pop]

#pxy <- p[, .(pi = weighted.mean(x = avg_pi, w = count_comparisons, na.rm=T)), by = pop]

#dxy <- d[, .(pop = "HH_vs_hh", pi = weighted.mean(x = avg_dxy, w = no_sites, na.rm=T))]
#xy <- rbind(pxy, dxy)

#xy
