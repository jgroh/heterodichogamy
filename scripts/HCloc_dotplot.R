library(data.table)

aln <- fread("~/workspace/heterodichogamy/HC_locus_structure/dotplots/Lakota1-6.2-7.1Mb_Pawnee-6.15-6.95Mb.megablast.csv")

PawneeTEs <- fread("~/workspace/heterodichogamy/HC_locus_structure/Pawnee_HCloc_TEanno.txt", col.names = c("chrom", "annotator", "type", "start", "stop", "score", "strand", "idk", "info"))
Lakota1TEs <- fread("~/workspace/heterodichogamy/HC_locus_structure/Lakota_v1_HCloc_TEanno.txt", col.names = c("chrom", "annotator", "type", "start", "stop", "score", "strand", "idk", "info"))

PawneeTEs[, hap := 'Pawnee (h)']
Lakota1TEs[, hap := "Lakota1 (H)"]
HC_TEs <- rbind(PawneeTEs, Lakota1TEs)

HC_TE_counts <- HC_TEs[, .N, by = type]
setkey(HC_TE_counts, "N")
HC_TEs[, type := factor(type, level = HC_TE_counts[,type])]


ggplot(HC_TEs, aes(x = type, fill = type)) + geom_bar(stat = 'count') + 
  facet_wrap(~hap) + 
  theme_classic() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_blank(),
        legend.position = c(0.7, 0.67),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.key.size = unit(.35, 'cm'),
        strip.text = element_text(size = 14)) + 
  labs(x = "", y = "Count", fill = "")



aln[, query_start := query_start + 6200000]
aln[, query_end := query_end + 6200000]

aln[, subject_start := subject_start + 6150000]
aln[, subject_end := subject_end + 6150000]

# make plot of broader region
ggplot(aln) + 
  geom_segment(aes(x = query_start, xend = query_end, 
                   y = subject_start, yend = subject_end),
               linewidth = 0.4)  + 
  geom_hline(yintercept = 6462000) + 
  geom_hline(yintercept = 6662000)  + 
  geom_vline(xintercept = 6500000) + geom_vline(xintercept = 7000000)
  
# make plot of just HC locus
ggplot(aln[subject_start > 6462000 & subject_end < 6662000 & query_start > 6500000 & query_start < 7000000]) + 
  geom_segment(aes(x = query_start, xend = query_end, 
                   y = subject_start, yend = subject_end),
               linewidth = 0.4)   +
  theme_classic() + 
  theme(aspect.ratio = 2/5) + 
  geom_rect(data = Lakota1TEs[V3 == 'Gypsy_LTR_retrotransposon'], aes(xmin = V4, xmax = V5, ymin=6462000, ymax=6662000), alpha = 0.2)
