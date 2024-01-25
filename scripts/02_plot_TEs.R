library(data.table)
library(ggplot2)
library(magrittr)
library(viridis)
library(cowplot)


cnames <- unlist(strsplit("Chr,Source,Type,Start,End,Score,Strand,Phase,Attributes", split = ','))

Ccat_TE <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_cathayensis/BNU/Carya_cathayensis_BNU_v1.fa.mod.EDTA.TEanno.gff3", skip = 'Chr', col.names = cnames)
Csin_TE <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_sinensis/Carya_sinensis_BNU_v1_chrs.fa.mod.EDTA.TEanno.gff3", skip = 'Chr', col.names = cnames)

CillLak1_TE <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_illinoinensis/Lakota_v1.fna.mod.EDTA.TEanno.gff3", skip = "CM0", col.names = cnames)
CillPaw1_TE <- fread("~/workspace/heterodichogamy/Carya_genome_assemblies/Carya_illinoinensis/Pawnee_RefSeq.fna.mod.EDTA.TEanno.gff3", skip = "NC_", col.names = cnames)

Ccat_TE[, id := 'cat']
Csin_TE[, id := 'sin']
CillLak1_TE[, id := 'lak1']
CillPaw1_TE[, id := 'paw']


TE <- rbindlist(list(Ccat_TE, Csin_TE, CillLak1_TE, CillPaw1_TE))



TE[, HCloc := 'Whole genome']

TE[(id == 'cat' & Chr == 'Chr04' & Start > 7186465 & End < 7419255) |
     (id == 'sin' & Chr == 'Chr04' & Start > 6074902 & End < 6421582) | 
     (id == 'lak1' & Chr == 'CM031828.1' & Start > 6507677 & End < 6952981) |
       (id == 'paw' & Chr == 'NC_056755.1' & Start > 6461560 & End < 6658636),
   HCloc := "G-locus"] 

TE[(id == 'cat' & Chr == 'Chr04' & ( (Start > 7186465 - 150000 & End < 7186465) | (Start > 7419255 & End < 7419255 + 150000) ) ) |
     (id == 'sin' & Chr == 'Chr04' & ( (Start > 6074902 - 150000 & End < 6074902) | (Start > 6421582 & End < 6421582 + 150000) ) ) |
     (id == 'lak1' & Chr == 'CM031828.1' & ( (Start > 6507677 - 150000 & End < 6507677) | (Start > 6952981 & End < 6952981 + 150000) ) ) |
     (id == 'paw' & Chr == 'NC_056755.1' & ( (Start > 6461560 - 150000 & End <  6461560) | ( Start > 6658636 & End < 6658636 + 150000) ) ),
   HCloc := "G-locus bordering regions"] 



# prepare data
TE[!id %in% c('cat', 'sin'), clade := 'pecan']
TE[id %in% c('cat', 'sin'), clade := 'EA clade']
TE[, clade := factor(clade, levels = c("pecan", "EA clade"))]



TE[id %in% c('cat','paw'), hap := 'g']
TE[id %in% c('sin','lak1'), hap := 'G']


# reorder factor
wg_counts <- TE[, .N, by = Type]
setkey(wg_counts, N)
TE[, Type := factor(Type, levels = wg_counts[, Type])]



pecan_main_TE_plt <- ggplot(TE[clade == 'pecan'], aes(x = Type, fill = Type)) + 
  geom_bar(stat = 'count') + 
  facet_grid(HCloc ~ hap, scales = 'free_y') + 
  theme_classic() + 
  labs(y = 'Count', fill = '', title = 'Pecan') +
  theme(aspect.ratio = 1, 
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

pecan_main_TE_plt 


leg <- ggplot(TE[HCloc == 'G-locus' & clade == 'pecan' & hap == 'G'], aes(x = Type, fill = Type)) + 
  geom_bar(stat = 'count') + 
  facet_grid(clade ~ hap) + 
  theme_classic() + 
  labs(x = 'Count', fill = '', title = 'Whole genome') +
  theme(aspect.ratio = 1, 
        #legend.position ='bottom',
        #legend.box = 'horizontal',
        legend.text = element_text(size = 10),
        legend.key.size = unit(.5, 'cm'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

leg

EA <- ggplot(TE[clade == 'EA clade' & HCloc == 'G-locus'], aes(x = Type, fill = Type)) + 
  geom_bar(stat = 'count') + 
  facet_grid(~ hap, scales = 'free') + 
  theme_classic() + 
  labs(y = 'Count', fill = '', title = 'East Asian clade') +
  theme(aspect.ratio = 1, 
        #legend.position ='bottom',
        #legend.box = 'horizontal',
        legend.text = element_text(size = 10),
        legend.key.size = unit(.5, 'cm'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank())

EA



