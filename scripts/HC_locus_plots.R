library(data.table)
library(readxl)
library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)
library(magrittr)

# 
meta <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")

# ---- read GWAS data -----
gbs <- fread("~/workspace/heterodichogamy/pecan/GBS_cultivars_gemma.assoc.txt")
gbs[, N:= seq_len(.N)]

wgs <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set_Lakota_v1.assoc.txt")
wgs[, N:= seq_len(.N)]

# ----- read pixy divergence data
# divergence mapped to Lakota, CDS

dxy <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/HH_vs_hh_Lakota1_CDS_dxy.txt")

# individual level heterozygosity.
# for heterozygotes this is just Dxy of haplotypes

ind_pi <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/ind_het_Lakota_v1_CDS_pi.txt")
ind_pi <- merge(ind_pi, meta[, .(ID, pop = variety, phenotype, genotype)], by = "pop")



# divergence mapped to Chandler, fixed windows
#dxy500 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/HH_vs_hh_500_dxy.txt")
#dxy5kb <- fread("~/workspace/heterodichogamy/pecan/WGS_data/results/pixy/HH_vs_hh_5000_dxy.txt")

#ggplot(dxy500[window_pos_1 > 6e6 & window_pos_2 <7e6],
#       aes(x = window_pos_1, y = avg_dxy)) + 
#  geom_point()




# -----  read phenotypes of WGS data
fam <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.fam")
fam <- fam[, .(V1, V6)]
setnames(fam, c("V1", "V6"), c("sample", "phenotype"))

#wgs_samples1 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples.tsv")
#setnames(wgs_samples1, 'run', 'sample')

#wgs_samples2 <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_reassigned.tsv")
#setnames(wgs_samples2, 'run', 'sample')

#geno <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples_geno.txt", header = F, col.names = c('run', 'geno'))

#my_samples <- as.data.table(read_excel("~/workspace/heterodichogamy/Manuscript/Supplementary_tables.xlsx", sheet = 'NovaSeq'))

#wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wg)


# ----- read LD data -----
ld <- fread("~/workspace/heterodichogamy/pecan/WGS_data/WGS_cultivars.geno.ld")
setnames(ld, c('chr', 'pos1', 'pos2', 'n_indv', 'r2'))

# ----- read HC locus genes
gff <- fread("~/workspace/heterodichogamy/HC_locus_structure/Lak1_mRNA_chr4_6.3-7.15Mb.txt")
#gn <- fread("~/workspace/heterodichogamy/HC_locus_structure/HC_locus_genes.txt")
#gff <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_genes_6.45-6.65Mb.GFF3")

# -------- library IDs of Carya spp

sra <- fread("~/workspace/heterodichogamy/Carya_spp/Carya_15sp_libraries.txt")

# --  read coverage of WGS data
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

#cvg <- merge(cvg, wgs_samples2)
#cvg <- merge(cvg, geno[, .(sample = run, geno)])

#cvg <- merge(cvg, my_samples[, .(sample = ID, phenotype, variety)])
#cvg[phenotype == 'protandrous', genotype := 'hh']
#cvg[phenotype == 'protogynous' & variety != 'Mahan', genotype := 'Hh']
#cvg[phenotype == 'protogynous' & variety == 'Mahan', genotype := 'HH']

cvg <- merge(cvg, fam)






# ------ read coverage of 15 sp ------
cvg_spp <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Carya_spp/results/coverage_Pawnee/", 
                                       pattern = '*.txt.gz', full.names = T),
                            function(x){
                              z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                              z[, sample := gsub(".txt.gz", "", basename(x))]
                              return(z)
                            }))
cvg_spp_nrm <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/Carya_spp/results/coverage_Pawnee/", 
                                       pattern = '*_norm.txt', full.names = T), 
                            function(x) {
                              z <- fread(x, col.names = 'avg_cvg')
                              z[, sample := gsub("_norm.txt", "", basename(x))]
                              return(z)
                            }
))
cvg_spp <- merge(cvg_spp, cvg_spp_nrm)


# ----- Plot GWAS

# # whole genome 
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


# ------ GWAS mapped to Lakota1 -------
max(wgs[, -log10(p_lrt)])

#CM031828.1:6500000-6700000
gff[, HC_loc := ifelse(V4 >= 6507677 & V5 <= 6952981, 'HC', 'BG')]

gwas_wg_Hloc <- ggplot(wgs[chr == 'CM031828.1' & ps > 6.4e6 & ps < 7.05e6]) +
  labs(y = expression(-log[10](P)),
       x = 'Chr 4 H hap (Mb) ') +
  scale_x_continuous(breaks = seq(6.4e6,7e6, length.out = 4), labels = c("6.4", '6.6', '6.8', '7.0')) +
  theme_classic() +
  theme(aspect.ratio = .3,
     #    plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2),
        #legend.position = c(0.15, 0.85),
        legend.position = 'none') +
  #geom_rect(data = gff, aes(xmin = V4, xmax = V5, ymin = 14, ymax = 16, fill = HC_loc)) + 
  geom_point(size =0.6, aes(x = ps, y = -log10(p_lrt))) +  

  geom_rect(data = gff[V4 > 6.4e6 & V5 < 7.05e6], aes(xmin = V4, xmax = V5, ymin = 3, ymax = 5, fill = HC_loc)) + 
  
  scale_fill_manual(values = c("gray", 'red')) 
 # geom_rect(data = gn, aes(xmin = Start, xmax = Stop, ymin=18, ymax=22), alpha = 0.3)

#geom_vline(aes(xintercept = 6461400)) + 
#geom_vline(aes(xintercept = 6663000)) +
gwas_wg_Hloc

#ggsave(filename = '~/workspace/heterodichogamy/Manuscript/figures/main/Jhindsii_pvals.pdf', plot = pval_plt, width = 7, height = 2, units = 'in')





# H-locus mapped to Pawnee
# max(wgs[, -log10(p_lrt)])
# #ps > 6.35e6 & ps < 6.77e6
# gwas_wg_Hloc <- ggplot(wgs[chr == 'CM031812.1' & ps > 6.45e6 & ps < 6.68e6]) +
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#         # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         #axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         legend.position = 'none') + 
#   geom_rect(data = gn, aes(xmin = Start, xmax = Stop, ymin=18, ymax=22), alpha = 0.3) 

  #geom_vline(aes(xintercept = 6461400)) + 
  #geom_vline(aes(xintercept = 6663000)) +

gwas_wg_Hloc

# ---- -Examine the other locations ----TEs???

# cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))
# CillPaw1_TE <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_illinoinensis/Pawnee_RefSeq.fna.mod.EDTA.TEanno.gff3", skip = "NC_", col.names = cnames)
# 
# ggplot(wgs[chr == 'CM031809.1'& ps > 2.107e7 & ps < 2.11e7 ]) +
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#         # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         #axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         legend.position = 'none')  + 
#   geom_rect(data = CillPaw1_TE[Chr == 'NC_056755.1' & Start > 2.107e7 & End < 2.11e7], aes(xmin=Start, xmax = End, ymin = 10, ymax = 15), fill = 'gray') + 
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) 
# 
# ggplot(wgs[chr == 'CM031811.1' & ps > 1.175e7 & ps < 1.178e7]) +
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#         # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         #axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         legend.position = 'none')  + 
#   geom_rect(data = CillPaw1_TE[Chr == 'NC_056754.1' & Start >1.175e7 & End < 1.178e7], aes(xmin=Start, xmax = End, ymin = 10, ymax = 15, fill= Type)) + 
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) 
# 
# ggplot(wgs[chr == 'CM031811.1' & ps > 2.08e7 & ps < 2.1e7 ]) +
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#         # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         #axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         legend.position = 'none')  + 
#   geom_rect(data = CillPaw1_TE[Chr == 'NC_056754.1' & Start >2.08e7 & End < 2.1e7], 
#             aes(xmin=Start, xmax = End, ymin = 10, ymax = 15, fill = Type)) + 
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) 
# 
# ggplot(wgs[chr == 'CM031811.1' & ps >4.567e7 & ps < 4.57e7]) +
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) +
#   labs(y = expression(-log[10](P)),
#        x = '') +
#   theme_classic() +
#   theme(aspect.ratio = .3,
#         # plot.margin = margin(30,30,30,30, "pt"),
#         text = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         #axis.text.x = element_blank(),
#         axis.title.x = element_text(size = 12, vjust = -2),
#         axis.title.y = element_text(size = 12, vjust = 2),
#         #legend.position = c(0.15, 0.85),
#         #legend.position = 'none'
#         )  + 
#   geom_rect(data = CillPaw1_TE[Chr == 'NC_056754.1' & Start >4.567e7 & End < 4.57e7], 
#             aes(xmin=Start, xmax = End, ymin = 10, ymax = 15, fill = Type)) + 
#   geom_point(size =1, aes(x = ps, y = -log10(p_lrt))) 
#   




# # whole genome 
# gwas_wg <- ggplot(gbs, aes(x = N, y = -log10(p_wald), color = chr)) + 
#          geom_point(size =1) + 
#          scale_color_manual(values = rep(c("tan", "purple4"), 8)) + 
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



# ----- zoomed in on GWAS peak -----
# gwas_peak <- ggplot(gbs[chr == 'CM031812.1' & ps > 6.4e6 & ps < 6.75e6], aes(x = ps, y = -log10(p_wald))) + 
#   geom_point() + 
#   #scale_color_manual(values = rep(c("black", "grey"), 8)) + 
#   theme_classic() +
#   labs(y = expression(-log[10](P)), 
#        x = '') +
#   theme(aspect.ratio = .3,
#       #plot.margin = margin(30,30,30,30, "pt"),
#       text = element_text(size = 12), 
#       axis.text = element_text(size = 12),
#       axis.text.x = element_blank(),
#       axis.title.x = element_text(size = 12, vjust = -2),
#       axis.title.y = element_text(size = 12, vjust = 2),
#       #legend.position = c(0.15, 0.85), 
#       legend.position = 'none')
#   #annotate("rect", xmin = 6534137, xmax = 6534830, ymin = 0, ymax = 25,
#   #          alpha = .8,fill = '#94435b') + 
#   #annotate("rect", xmin = 6632000, xmax = 6634000, ymin = 0, ymax = 25,
#   #         alpha = .8,fill = '#94435b') 




# ---- calculate coverage in  windows
#Pawnee
#st <- 6000000
#en <- 7000000

st <- 6400000
en <- 7050000

# --- 1kb
cvg[, window := cut(position, breaks = seq(st, en+1000, by = 1000), labels = seq(st, en, by = 1000), include.lowest =T), by = sample]

cvg1kb <- cvg[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg, phenotype)]
cvg1kb[, window := as.numeric(as.character((window)))]

cvg1kb[, nrm_cvg := coverage/avg_cvg]

# looks like there presence/absence variation
# ggplot(cvg1kb[window >6400000 & window < 6700000], 
#        aes(x = window, y = nrm_cvg, group = sample, color = as.factor(phenotype))) + geom_line() +
#   geom_line(data = cvg1kb[sample == 'SRR15911540' & window >6400000 & window < 6700000]) 


# --- 10bp
# cvg[, window10 := cut(position, breaks = seq(st, en+10, by = 10), labels = seq(st, en, by = 10), include.lowest =T), by = sample]
# 
# cvg10 <- cvg[, .(coverage = mean(coverage)), by = .(sample,window10, avg_cvg, type, phenotype, genotype)]
# cvg10[, window := as.numeric(as.character((window10)))]
# 
# cvg10[, nrm_cvg := coverage/avg_cvg]




# ---- -calculate coverage for 15sp
# cvg_spp[, window := cut(position, breaks = seq(st, en+1000, by = 1000), labels = seq(st, en, by = 1000), include.lowest =T), by = sample]
# 
# cvg_spp1kb <- cvg_spp[, .(coverage = mean(coverage)), by = .(sample,window, avg_cvg)]
# cvg_spp1kb[, window := as.numeric(as.character((window)))]
# 
# cvg_spp1kb[, nrm_cvg := coverage/avg_cvg]
# cvg_spp1kb <- merge(cvg_spp1kb, sra[, .(sample = Run, Species, putative_geno)])

# ggplot(cvg_spp1kb[window > 6.42e6 & window < 6.67e6], aes(x = window, y = nrm_cvg, group = sample, color = putative_geno))  + 
#   geom_line() + 
#   theme_classic() +
#   theme(
#     aspect.ratio = .2,
#     #plot.margin = margin(30,30,30,30, "pt"),
#     text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     #axis.text.x = element_blank(),
#     axis.title.x = element_text(size = 12),
#     axis.title.y = element_text(size = 12, vjust = 2),
#     legend.title = element_blank(),
#     legend.position = c(0.9, 0.95),
#     legend.margin=margin(c(0,0,0,0)),
#     legend.key.size = unit(0.3, 'cm'),
#     #plot.margin = unit(c(0,0,0,2), 'cm')) + # if bottom
#     plot.margin = unit(c(0,0,-4,2), 'cm'))

# ====== Commented out sections below were used to identify the SV genotypes with publicly available Carya reseq
# for which variety labels appeared to be mismatched
# but superseded

# ----- color by phenotype *** note that there is clear evidence that some phenotypes were mismatched in the wgs data set ****
#ggplot(cvg1kb[window >6640000 & window < 6642500], 
#       aes(x = window, y = nrm_cvg, 
#           group = sample, color = as.factor(phenotype))) + geom_line()

# which samples don't match???
# cvg1kb[window >6640000 & window < 6642500 & nrm_cvg < 0.8 & nrm_cvg > 0.2, unique(type)]

# ----- color by phenotype *** FOR SUBSET OF FILTERED INDIVIDUALS BASED ON RELATEDNESS TO GBS SAMPLES ****
# Note this section was prev

# type %in% c("Barton", "Dependable", "Excel", "Farley","Kiowa", "Mahan", "Oconee", "Pawnee", "Shoshoni","Sioux", "Sumner", "VC1-68"
#ggplot(cvg1kb[window >6600000 & window < 6642500], 
#       aes(x = window, y = nrm_cvg, 
#           group = sample, color = as.factor(phenotype))) + geom_line()

#highconf_samples <- c("Barton", "Dependable", "Excel", "Farley",
#                      "Kiowa", "Mahan", "Oconee", "Pawnee", "Shoshoni",
#                      "Sioux", "Sumner", "VC1-68")
#cvg1kb[!type %in% highconf_samples, genotype := "NA"]
#cvg1kb[type == 'Mahan', genotype := 'HH']
#cvg1kb[type != 'Mahan' & phenotype == 'protogynous', genotype := 'Hh']
#cvg1kb[phenotype == 'protandrous', genotype := 'hh']

# ----- assign SV SV genotype *** for all samples, including those which phenotypes are off for 
# ggplot(cvg1kb[window >6640000 & window < 6642500], 
#        aes(x = window, y = nrm_cvg, 
#            group = sample)) + geom_line()
# SV_HH <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg < 0.1, unique(sample)]
# cvg1kb[sample %in% SV_HH, SV := 2]
# 
# SV_Hh <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg > 0.1 & nrm_cvg < 0.8, unique(sample)]
# cvg1kb[sample %in% SV_Hh, SV := 1]
# 
# SV_hh <- cvg1kb[window >6640000 & window < 6642500 & nrm_cvg > 0.8, unique(sample)]
# cvg1kb[sample %in% SV_hh, SV := 0]
# 
# # check that the SV is scored correctly
# ggplot(cvg1kb[window >6640000 & window < 6642500], 
#        aes(x = window, y = nrm_cvg, 
#            group = sample, color = as.factor(SV))) + geom_line()


#cvg1kb[window >6400000 & window < 6750000]
# # now zoom out again
# plot3 <- ggplot(cvg1kb[window >6520000 & window < 6700000 & SV == 0 & sample == 'SRR15911540'],
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
# plot3
# 
# # ----- write out genotypes -----
# cvg1kb[SV==0, genotype := 'hh']
# cvg1kb[SV==1, genotype := 'Hh']
# cvg1kb[SV==2, genotype := 'HH']
# 
# fwrite(cvg1kb[ window == min(window), .(sample, genotype)], file = 'wgs_samples_geno.txt',
#        col.names = F, row.names = F, quote = F, sep = "\t")
# 


# ========= MAIN PLOTS =======

# dxy500[, avgDxy := mean(avg_dxy, na.rm=T)]
# dxy500[, qnt0.99 := quantile(avg_dxy, probs = c(0.99), na.rm=T)]
# dxy500[, outlier := ifelse(avg_dxy > qnt0.99, 1, 0)]

# dxy[, avgDxy := mean(avg_dxy, na.rm=T)]
# dxy[, qnt0.99 := quantile(avg_dxy, probs = c(0.99), na.rm=T)]
# dxy[, outlier := ifelse(avg_dxy > qnt0.99, 1, 0)]


# ----- dxy -------

# dxy_plt <- ggplot(dxy[chromosome == 'CM031828.1' & window_pos_1 > 5e6 & window_pos_2 < 8e6], 
#                    aes(x = (window_pos_1 + window_pos_2)/2, y = avg_dxy, color = as.factor(outlier))) + 
#   geom_rect(data = dxy[chromosome == 'CM031828.1' & window_pos_1 > 5e6 & window_pos_2 < 8e6],
#             aes(xmin=-Inf, xmax=Inf, ymin=0, ymax = qnt0.99), fill = 'gray95', 
#             alpha = 0.5, color = 'NA') + 
#   geom_point(size =1) +  
#   geom_hline(data = dxy[chromosome == 'CM031828.1' & window_pos_1 > 5e6 & window_pos_2 < 8e6], 
#              aes(yintercept = avgDxy), color = 'darkgray', linetype = 2) +
#   scale_color_manual(values = c("black", 'red')) + 
#   labs(y = expression('D'[xy]), x = 'Pawnee chr 4 (Mb)') + # if bottom
#   #labs(x = "", y = expression('D'[xy]), x = x = '') + # if top
#   
#   
#   #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = rep("",5))  + # if top
#   scale_x_continuous(breaks = seq(6.4e6, 7e6, length.out = 3), labels = seq(6.4, 7, length.out = 3))  + # if bottom
#   
#   #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels =  NULL)  +
#   theme_classic() + 
#   #geom_vline(xintercept = 6462000, color ='red') +
#   #geom_vline(xintercept = 6662000, color = 'red') +
#   theme(aspect.ratio = 0.2, 
#         legend.position = 'none', 
#         legend.key.size = unit(0.5, 'cm'),
#         axis.text = element_text(size = 10), 
#         axis.title = element_text(size = 12), 
#         #plot.margin = unit(c(0,0,-4,2), 'cm')) # if top
#         plot.margin = unit(c(0,0,0,2), 'cm'))  # if bottom

# plot using individual level heterozygosity 

dxy_data <- rbind(dxy[, .(pop = 'HH_vs_hh', chromosome, window_pos_1, window_pos_2, avg_dxy)],
                  ind_pi[genotype == 'Hh',  .(pop, chromosome, window_pos_1, window_pos_2, avg_dxy = avg_pi)])

dxy_avg_gw <- dxy_data[, .(dxy = mean(avg_dxy, na.rm=T)), by = .(chromosome, window_pos_1, window_pos_2)] 

# mean and outliers
dxy_avg_gw[, avgDxy := weighted.mean(dxy,w = window_pos_2 - window_pos_1, na.rm=T)]
dxy_avg_gw[, qnt0.95 := quantile(dxy, probs = c(0.95), na.rm=T)]
dxy_avg_gw[, outlier := ifelse(avg_dxy > qnt0.95, 1, 0)]


dxy_plot_data <- dxy_avg_gw[chromosome == 'CM031828.1' & window_pos_1 > 6.4e6 & window_pos_2 < 7.05e6] 

  
dxy_plt <- ggplot(dxy_plot_data,
       aes(x = (window_pos_1 + window_pos_2)/2, y = dxy)) + 
  geom_rect(aes(xmin=-Inf, xmax= Inf, ymin=0, ymax = qnt0.95), fill = 'gray95', alpha = 0.5, color = 'NA') + 
  geom_hline(aes(yintercept = avgDxy), linetype = 2, color = 'darkgray') +
  geom_point(size =1) +  
  scale_x_continuous(breaks = seq(6.4e6, 7e6, length.out = 3), labels = seq(6.4, 7, length.out = 3))  + 
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels =  NULL)  +
  theme_classic() + 
  labs(y = expression('D'[xy]), x = 'Chr 4 H hap (Mb)') + # if bottom
  theme(aspect.ratio = 0.2, 
        legend.position = 'none', 
        legend.key.size = unit(0.5, 'cm'),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        #plot.margin = unit(c(0,0,-4,2), 'cm')) # if top
        plot.margin = unit(c(0,0,0,2), 'cm')) # if bottom
        
dxy_plt






ld_plt <- ggplot(ld[pos1 > 6.42e6 & pos1 < 6.7e6 & pos2 >  6.42e6 & pos2 < 6.7e6], 
                 aes(x = pos2, y = pos1, fill = r2)) + 
  geom_tile(width = 200, height = 200) +
  #annotate("rect", xmin = 31868751, xmax = 31870103, 
  #         ymin = 31868751, ymax = 31870103,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  #annotate("rect", xmin = 31884268, xmax = 31887072, 
  #         ymin = 31884268, ymax = 31887072,
  #         alpha = .9,fill = 'NA', color = 'black') + 
  scale_fill_viridis() + 
  labs(x = "Pawnee chr 4 (Mb)", y = "Pawnee chr 4 (Mb)", fill = expression(r^2)) +
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = NULL )  +
  scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = sprintf("%.2f", seq(6.45, 6.65, length.out = 5))) +
  
  scale_y_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = sprintf("%.2f", seq(6.45, 6.65, length.out = 5))) +
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
  annotate("rect", xmin = 6461560, xmax = 6658636, 
           ymin = 6461560, ymax = 6658636,
           alpha = .9,fill = 'NA', color = 'black') 

#6461560 & End < 6658636

ld_plt

#window >6420000 & window < 6700000
#cvg1kb[phenotype == 0, genotype := 'hh']
#cvg1kb[phenotype == 1, genotype := 'H?']
#cvg1kb[sample == 'SRR15911533', genotype := 'HH']
#cvg1kb[sample == 'CILL_LW1_11', genotype := 'HH']
#cvg1kb[, genotype := factor(genotype, levels = c("HH", 'H?', 'hh'))]
cvg_plt <- ggplot(cvg1kb[window >6.4e6 & window < 7.05e6],
            ) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample,color = genotype), linewidth = 0.5) +
  #geom_line(data = cvg1kb[variety == 'Apachee' & window >6420000 & window < 6700000], aes(x = window, y = nrm_cvg), color = 'black') +
  
  #scale_y_continuous(limits = c(0, 3)) +
  scale_color_manual(values = c('turquoise4', 'maroon', 'tan')) +
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = rep("",5))  + # if top
  
  #scale_x_continuous(breaks = seq(6.45e6, 6.65e6, length.out = 5), labels = seq(6.45, 6.65, length.out = 5))  + # if bottom
  scale_y_continuous(limits = c(0,3)) +
  #geom_vline(aes(xintercept = 6.662e6)) +
  theme_classic() +
  #labs(y = 'Normalized\nread depth', x = 'Pawnee chr 4 (Mb)', color = '') + # if bottom
  labs(y = 'Normalized\nread depth', x = 'Chr 4 H hap (Mb)', color = '') + # if top
  scale_x_continuous(breaks = seq(6.4e6,7e6, length.out = 4), labels = c("6.4", '6.6', '6.8', '7.0')) +
  
  theme(aspect.ratio = .3,
        #plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2),
        #legend.position = c(0.15, 0.85),
        legend.position = c(0.9, 0.8),
        legend.key.size = unit(0.3, 'cm'),
        
  ) #+
  #geom_rect(data = gff, aes(xmin = V4, xmax = V5, ymin = 14, ymax = 16, fill = HC_loc)) + 
  

cvg_plt

gg <- plot_grid(gwas_wg_Hloc, cvg_plt, dxy_plt, ncol=1, align='v')
gg


#grid.arrange(ld_plt, cvg_plt, ncol = 1)
plot_grid(dxy_plt, ld_plt, cvg_plt, ncol=1, align='v')
plot_grid(cvg_plt, ld_plt, dxy_plt, ncol=1, align='v')




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





# ----- coverage plots for each gene individually -----

# PPI2
ggplot(cvg10[window > 6461404 -100 & window < 6464894 + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'protein phosphatase inhibitor 2-like') + # if top
  theme(
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122306596', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122306596', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))


# uncharacterized, LOC122306891

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122306891", Start] - 100 & window < gn[`Gene symbol` == "LOC122306891", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122306891 (uncharacterized)') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122306891', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122306891', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# EMS1, LOC122308212

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308212", Start] - 100 & window < gn[`Gene symbol` == "LOC122308212", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'leucine-rich repeat receptor protein kinase EMS1-like') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308212', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308212', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))




# CEN, LOC122308000

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308000", Start] - 100 & window < gn[`Gene symbol` == "LOC122308000", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'CEN-like protein 1') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308000', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308000', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# XBAT33, LOC122307679

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307679", Start] - 100 & window < gn[`Gene symbol` == "LOC122307679", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'E3 ubiquitin-protein ligase XBAT33-like') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
   # plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307679', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307679', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))




# RPD1, LOC122308480

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308480", Start] - 100 & window < gn[`Gene symbol` == "LOC122308480", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'protein ROOT PRIMORDIUM DEFECTIVE 1') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308480', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308480', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# RHOMBOID2

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308345", Start] - 100 & window < gn[`Gene symbol` == "LOC122308345", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'RHOMBOID-like protein 2') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308345', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308345', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# uncharacterized, LOC122306943

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122306943", Start] - 100 & window < gn[`Gene symbol` == "LOC122306943", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122306943 (uncharacterized)') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122306943', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122306943', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))


# putative non-specific lipid-transfer protein 14, LOC122306944

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122306944", Start] - 100 & window < gn[`Gene symbol` == "LOC122306944", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'non-specific lipid-transfer protein 14') + # if top
  theme(
    aspect.ratio = .2,
    legend.position = 'None',
    plot.margin = margin(r = 20, l = 20, unit = 'pt'),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    #plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122306944', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122306944', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))





# IQ-DOMAIN 23-like, LOC122307664

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307664", Start] - 100 & window < gn[`Gene symbol` == "LOC122307664", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'IQ-DOMAIN 23-like') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307664', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307664', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# FIL1-like, LOC122308254

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308254", Start] - 100 & window < gn[`Gene symbol` == "LOC122308254", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'FIL1') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308254', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308254', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))


# PMT16, LOC122307000

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307000", Start] - 100 & window < gn[`Gene symbol` == "LOC122307000", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'PMT16') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307000', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307000', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# DNA-directed RNA polymerases..., LOC122307855

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307855", Start] - 100 & window < gn[`Gene symbol` == "LOC122307855", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122307855') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307855', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307855', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# F-box protein PP2-A12-like, LOC122306681

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122306681", Start] - 100 & window < gn[`Gene symbol` == "LOC122306681", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'PP2-A12') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122306681', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122306681', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# F-box protein At4g00755-like, LOC122307040

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307040", Start] - 100 & window < gn[`Gene symbol` == "LOC122307040", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'At4g00755') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307040', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307040', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))





# SLK2, LOC122307039

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307039", Start] - 100 & window < gn[`Gene symbol` == "LOC122307039", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'SLK2') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307039', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307039', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))


#protein COFACTOR ASSEMBLY OF COMPLEX C SUBUNIT B CCB2, chloroplastic-like
# LOC122307672

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307672", Start] - 100 & window < gn[`Gene symbol` == "LOC122307672", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122307672') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307672', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307672', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# CRR1, LOC122307671

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307671", Start] - 100 & window < gn[`Gene symbol` == "LOC122307671", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122307671') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307671', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307671', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))


# BAG1, LOC122308043


ggplot(cvg10[window > gn[`Gene symbol` == "LOC122308043", Start] - 100 & window < gn[`Gene symbol` == "LOC122308043", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122308043') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122308043', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122308043', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))



# protein DETOXIFICATION 49-like, LOC122307012

ggplot(cvg10[window > gn[`Gene symbol` == "LOC122307012", Start] - 100 & window < gn[`Gene symbol` == "LOC122307012", Stop] + 100],
) + 
  geom_line(aes(x = window, y = nrm_cvg, group = sample, color = geno)) +
  scale_color_manual(values = c('maroon', 'tan', 'turquoise')) +
  
  scale_y_continuous(limits = c(0,3)) +
  theme_classic() +
  labs(y = 'Normalized\nread depth', x = '', color = '', title = 'LOC122307012') + # if top
  theme(
    aspect.ratio = .2,
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12, vjust = 2),
    legend.title = element_blank(),
    plot.title = element_text(face = 'italic'),
    legend.key.size = unit(0.3, 'cm')) +
  geom_segment(data = gff[V3=='mRNA' & grepl('LOC122307012', V9)], aes(x=V4, xend=V5, y = 3, yend = 3)) +
  geom_rect(data = gff[V3 == 'CDS' & grepl('LOC122307012', V9)], aes(xmin = V4, xmax = V5, ymin=2.5, ymax = 3))




