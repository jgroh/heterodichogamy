library(data.table)
library(readxl)
library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)
library(magrittr)

# ===== Read Data =====
meta <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")

# ---- gwas ----
gbs <- fread("~/workspace/heterodichogamy/pecan/GBS_cultivars_gemma.assoc.txt")
gbs[, N:= seq_len(.N)]

wgs <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_Lakota_v1.assoc.txt")
wgs[, N:= seq_len(.N)]

# ----- divergence (in CDS, mapped to Lakota v1) -----
# ----- homozygotes
#dxy <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/HH_vs_hh_Lakota1_CDS_dxy.txt")
#dxy <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/GG_vs_gg_dxy.txt")
dxy <- fread("~/workspace/heterodichogamy/01_GC_locus_structure/Lakota_v1_Pawnee_CDS_dxy.txt")

# individual level heterozygosity.
# for heterozygotes at H locus this is Dxy between the haplotypes

ind_pi <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/ind_het_Lakota_v1_CDS_pi.txt")
ind_pi <- merge(ind_pi, meta[, .(ID, pop = variety, phenotype, genotype)], by = "pop")


# ----- read heterozygosity of other Carya spp
spp <- fread("~/workspace/heterodichogamy/Carya_spp/results/pixy/Carya_spp_pi_at_pecan_SNPS_pi.txt")
spp_het <- spp[, .(p = mean(avg_pi, na.rm=T)), by = pop]

# calculate standard errors
spp_N <- spp[no_sites == 1, .N, by = pop]
spp_het <- merge(spp_het, spp_N)

ENA <- c("aquatica", "cordiformis", "floridana", "glara", "laciniosa", 
         "myristiciformis", "ovata", "palmeri", "texana", "tomentosa")
EA <- c("cathayensis", "dabieshanensis", "hunanensis", 'ovata_1', 'ovata_2', 'kweichowensis', 'poilanei', "tonkinensis_1", "tonkinensis_2")

spp_het[pop %in% ENA, clade := 'North American']
spp_het[pop %in% EA, clade := 'East Asian']
spp_het[pop == 'glara', pop := 'glabra']
spp_het[pop == 'tonkinensis_2', pop := 'tonkinensis']
setkey(spp_het, clade)
lv2 <- rev(unlist(spp_het[, pop]))
spp_het[,  pop := factor(pop, levels = lv2)]

# to add rug 
het_ids_NA <- spp_het[p > 0.2 & clade == 'North American', pop]
het_ids_EA <- spp_het[p > 0.2 & clade == 'East Asian', pop]



# ----- read LD data -----
#ld <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.geno.ld")
ld <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/LD_Lakota_v1.geno.ld")

setnames(ld, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))

# read TEs
cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
Lak1_TEs <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_illinoinensis/Lakota_v1.fna.mod.EDTA.TEanno.gff3", skip = "CM0", col.names = cnames)


# ----- read HC locus genes -----
gff <- fread("~/workspace/heterodichogamy/HC_locus_structure/Lak1_mRNA_chr4_6.3-7.15Mb.txt")

# ------  read coverage of WGS data -----
wgs_cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_v1/", 
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))

cvg_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/pecan/WGS_data/results/coverage_Lakota_v1///", 
                                       pattern = '*_norm.txt', full.names = T), 
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }
))

cvg <- merge(wgs_cvg, cvg_nrm)
cvg <- merge(cvg, meta[, .(sample = ID, variety, phenotype, genotype)])



# ===== Calculate windowed coverage =====

st <- 6400000
en <- 7050000

# --- 1kb
cvg[, window := cut(position, breaks = seq(st, en+1000, by = 1000), labels = seq(st, en, by = 1000), include.lowest =T), by = sample]

cvg1kb <- cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg, phenotype, genotype)]
cvg1kb[, window := as.numeric(as.character((window)))]

cvg1kb[, nrm_cvg := coverage/avg_cvg]

# --- 10bp
cvg[, window := cut(position, breaks = seq(st, en+10, by = 10), labels = seq(st, en, by = 10), include.lowest =T), by = sample]

cvg10 <- cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg, phenotype, genotype)]
cvg10[, window := as.numeric(as.character((window)))]

cvg10[, nrm_cvg := coverage/avg_cvg]


# ===== Plots =====

# ----- GWAS ------

# ----- whole genome 
# gwas_wg <- ggplot(wgs[seq(1, .N, by = 1)], aes(x = N, y = -log10(p_lrt), color = chr)) +
#          geom_point(size = 0.8) +
#          scale_color_manual(values = rep(c("gray", "black"), 8)) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#        # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         legend.position = 'none')
# 
# gwas_wg
# ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/pecan_GWAS.pdf', plot = gwas_wg, width = 7, height = 2, units = 'in')


# ----- HC locus

# define variable to indicate whether gene is in locus or outside
gff[, HC_loc := ifelse(V4 >= 6507677 & V5 <= 6952981, 'HC', 'BG')]


gwas_Hloc <- ggplot(wgs[chr == 'CM031828.1' & ps > 6.4e6 & ps < 7.05e6]) +
  labs(y = expression(-log[10](P)),
       x = '') +
  scale_y_continuous(limits = c(-2,20)) +
  scale_x_continuous(breaks = seq(6.4e6,7e6, length.out = 4), labels = c("6.4", '6.6', '6.8', '7.0')) +
  theme_classic() +
  theme(aspect.ratio = .2,
        # plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none') +
  geom_point(size =0.6, aes(x = ps, y = -log10(p_lrt))) +  
  geom_rect(data = gff[V4 > 6.4e6 & V5 < 7.05e6], aes(xmin = V4, xmax = V5, ymin = -2, ymax = 0, fill = HC_loc)) + 
  scale_fill_manual(values = c("gray", 'red')) 

gwas_Hloc


# ----- Coverage -----
cvg1kb[genotype == 'HH', genotype := 'GG']
cvg1kb[genotype == 'Hh', genotype := 'Gg']
cvg1kb[genotype == 'hh', genotype := 'gg']

cvg10[genotype == 'HH', genotype := 'GG']
cvg10[genotype == 'Hh', genotype := 'Gg']
cvg10[genotype == 'hh', genotype := 'gg']


cvg_plt <- ggplot(cvg1kb[window >6.4e6 & window < 7.05e6],
#ggplot(cvg10[window >6727154 & window < 6728426],

) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample,color = genotype), linewidth = 0.5) +

  scale_color_manual(values = c('tan', 'maroon', 'turquoise4')) +

  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '') + # if top
  scale_x_continuous(breaks = seq(6.4e6,7e6, length.out = 4), labels = c("6.4", '6.6', '6.8', '7.0')) +
  
  theme(aspect.ratio = .2,
        #plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        #legend.position = c(0.15, 0.85),
        legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.35, 'cm'),
        
  ) 

cvg_plt

# ---- Examine coverage over possible Rho-GTPase ----
# rho gtpase window > 6754457 & window < 6785572
# uncharacterized: 6827960-6827454. 
transloc_cvg <- ggplot(cvg10[window > 6827454 & window < 6827960],
                  #ggplot(cvg10[window >6727154 & window < 6728426],
                  
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample,color = genotype), linewidth = 0.5) +
  
  scale_color_manual(values = c('tan', 'maroon', 'turquoise4')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '') + # if top
  #scale_x_continuous(breaks = seq(6.4e6,7e6, length.out = 4), labels = c("6.4", '6.6', '6.8', '7.0')) +
  
  theme(aspect.ratio = .2,
        #plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        #legend.position = c(0.15, 0.85),
        legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.35, 'cm'),
        
  )  + 
  geom_segment(data = Lak1_TEs[ Start > 6827454 & End < 6827960], 
               aes(x = Start, xend = End, y = 2, yend = 2, color= 'red')) + 
  geom_segment(data = Lak1_TEs[ Start > 6827454 & End < 6827960], 
               aes(x = Start, xend = End, y = 2, yend = 2)) 
  

transloc_cvg





# ----- Dxy -----

## 6507677 & End < 6952981
#ggplot(dxy[chromosome == 'CM031828.1'], aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) + 
#  geom_point()
#ggplot(dxy[chromosome == 'CM031828.1' & window_pos_1 > 6.4e6 & window_pos_2 < 7.05e6], aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) + 
#  geom_point()
dxy[chromosome == 'CM031828.1' & window_pos_1 > 6507677 & window_pos_2 < 6952981, weighted.mean(avg_dxy, w = window_pos_2 - window_pos_1, na.rm = T)]



#dxy_data <- rbind(dxy[, .(pop = 'HH_vs_hh', chromosome, window_pos_1, window_pos_2, avg_dxy)],
#                  ind_pi[genotype == 'Hh',  .(pop, chromosome, window_pos_1, window_pos_2, avg_dxy = avg_pi)])

#dxy_avg_gw <- dxy[, .(dxy = mean(avg_dxy, na.rm=T)), by = .(chromosome, window_pos_1, window_pos_2)] 

# mean and outliers
#dxy_avg_gw[, avgDxy := weighted.mean(dxy,w = window_pos_2 - window_pos_1, na.rm=T)]
#dxy_avg_gw[, qnt0.95 := quantile(dxy, probs = c(0.95), na.rm=T)]
#dxy_avg_gw[, outlier := ifelse(dxy > qnt0.95, 1, 0)]

dxy[, avgDxy := weighted.mean(avg_dxy, w = window_pos_2 - window_pos_1, na.rm=T)]
dxy[, thresh := quantile(avg_dxy, probs = c(0.99), na.rm=T)]
dxy[, outlier := ifelse(avg_dxy > thresh, 1, 0)]

dxy_plot_data <- dxy[chromosome == 'CM031828.1' & window_pos_1 > 6.4e6 & window_pos_2 < 7.05e6] 


dxy_plt <- ggplot(dxy_plot_data,
                  aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy)) + 
  geom_rect(aes(xmin=-Inf, xmax= Inf, ymin=0, ymax = thresh), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_hline(aes(yintercept = avgDxy), linetype = 2, color = 'darkgray') +
  geom_point(size = 0.8) +  
  scale_x_continuous(breaks = seq(6.4e6, 7e6, length.out = 4), labels = sprintf("%.1f", seq(6.4, 7, length.out = 4)))  + 
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels =  NULL)  +
  theme_classic() + 
  #scale_y_continuous(limits = c(-.05, 0.6)) + 
  scale_y_continuous(limits = c(-.004, 0.06)) + 
  labs(y = expression('D'[xy]), x = expression('Chr 4' ~ italic(G) ~ 'hap (Mb)')) + # if bottom
  geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_NA], aes(x = window_pos_1, y = NULL), color = 'blue', length = unit(.3, 'cm')) + 
  geom_rug(data = spp[avg_pi == 1 & pop %in% het_ids_EA], aes(x = window_pos_1, y = NULL), color = 'red', length = unit(.3, 'cm')) + 
  #geom_rug(data = spp[avg_pi == 0 & pop %in% het_ids_EA], aes(x = window_pos_1, y = NULL), color = 'green', length = unit(.3, 'cm')) + 
  
  theme(aspect.ratio = .2,
        #plot.margin = margin(0,0,0,0, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        #legend.position = c(0.15, 0.85),
        legend.position = c(0.9, 0.8),
        legend.key.size = unit(0.3, 'cm'),
  ) 


dxy_plt

plot_grid(gwas_Hloc, cvg_plt, dxy_plt, ncol=1, align='v')

gg <- plot_grid(gwas_Hloc, cvg_plt, dxy_plt, ncol=1, align='v')
gg

spp[pop %in% het_ids_NA & window_pos_1 < 6600000 & pop %in% het_ids_NA & !is.na(avg_pi) ]

spp[window_pos_1 < 6600000 & pop %in% het_ids_EA & !is.na(avg_pi),]




# ----- TE plot -----
TEs <- Lak1_TEs[Chr == 'CM031828.1' & Start > 6.4e6 & End < 7.05e6, .(Type, start = Start, end = End)]


# ----- average TEs -----
TE_cvg <- fread("~/workspace/heterodichogamy/HC_locus_structure/Lak1_Hloc_TE_cvg_1kb.txt", 
                col.names = c("chr", 'start', 'end', 'n', 'nbase', 'winsize', 'TEcvg'))
TE_cvg[, pos := end - 500]

#ggplot(TE_cvg, aes(x = start, y = TEcvg)) + geom_point() + geom_smooth(method = 'loess')
# 
# w = 10e3
# windows <- data.table(chr='CM031828.1', pos = seq(6.4e6, 7e6, by = w))
# windows[, start := pos - w/2][, end := pos + w/2]
# fwrite(windows[, .(chr, start, end)], file = '~/workspace/heterodichogamy/HC_locus_structure/Hloc_1kb_windows.bed', quote = F, row.names = F, col.names = F, sep = "\t")
# # 
# 
# # find all fragments that overlap grid windows
# setkey(TEs, start, end)
# 
# setkey(windows, start, end)
# overlaps <- foverlaps(windows, TEs, by.x=c("start","end"), by.y=c("start","end"))
# 
# # percent coverage of window
# overlaps[, TE_cvg := (pmin(end, i.end) - pmax(start, i.start))/w , by=.(start,end)]
# 
# # windows with no overlap get weight zero
# overlaps[is.na(TE_cvg), TE_cvg := 0]
# 
# # count TEs per window
# nTEs <- overlaps[, .N, by = .(pos)]

#ggplot(nTEs, aes(x =pos, y = N )) + geom_point() + geom_smooth(method = 'loess')


# # take max coverage across fragments in each window
# TE_plot_data <- overlaps[, .(TE_cvg = max(TE_cvg)), by = .(i.start, pos, i.end)]
# setnames(TE_plot_data , c("i.start", "i.end"), c("start", "end"))


# cvg_diff_HH_v_hh <- dcast(cvg1kb[genotype %in% c("hh", 'HH'), .(cvg = mean(nrm_cvg)), by= .(window, genotype)], window ~ genotype)
# cvg_diff_HH_v_hh[, cvg_diff := HH-hh]
# cvg_TEs <- merge(cvg_diff_HH_v_hh, TE_cvg[,.(window = pos, TEcvg)])
# 
# cvg_TEs
# 
# ggplot(cvg_TEs, aes(x = TEcvg, y = cvg_diff)) + geom_point()  + geom_smooth(method = 'lm') + theme_classic()
# 




#Start > 6507677 & End < 6952981
ld
ld_plt <- ggplot(ld, 
                 aes(x = pos1, y = pos2, fill = r2)) + 
  geom_tile(width = 100, height = 100) +
  #annotate("rect", xmin = 31868751, xmax = 31870103, 
  #         ymin = 31868751, ymax = 31870103,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  #annotate("rect", xmin = 31884268, xmax = 31887072, 
  #         ymin = 31884268, ymax = 31887072,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  #labs(x = "Pawnee chr 4 (Mb)", y = "Pawnee chr 4 (Mb)", fill = expression(r^2)) +
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = NULL )  +
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = sprintf("%.2f", seq(6.45, 6.65, length.out = 5))) +
  
  #scale_y_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = sprintf("%.2f", seq(6.45, 6.65, length.out = 5))) +
  theme_classic() + 
  #geom_hline(yintercept = 6462000, color ='red') +
  #geom_hline(yintercept = 6662000, color = 'red') +
  #geom_vline(xintercept = 6462000, color ='red') +
  #geom_vline(xintercept = 6662000, color = 'red') +
  theme(aspect.ratio = 1, 
        legend.position = c(0.1, 0.7), 
        legend.key.size = unit(0.3, 'cm'),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        #plot.margin = unit(c(0,-2,-2,2), 'cm')
  ) + 
  annotate("rect", xmin = 6507677, xmax = 6952981, 
           ymin = 6507677, ymax = 6952981,
           alpha = .9,fill = 'NA', color = 'black') 

#6461560 & End < 6658636

ld_plt


