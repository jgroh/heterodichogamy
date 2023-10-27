library(lme4)
library(data.table)

p <- fread("~/workspace/heterodichogamy/data/heterodichogamy_phenotypes.txt")
d <- fread("~/workspace/heterodichogamy/data/New_Stuke_catkin_buds_20Jul2023.tsv")
library(ggplot2)

cat_cols <- paste0("catkin", 1:16)
wal_cols <- paste0("walnut", 1:6)

d[, catkin_mean := rowMeans(.SD, na.rm=T), .SDcols = cat_cols]
d[, walnut_mean := rowMeans(.SD, na.rm=T), .SDcols = wal_cols]

mns <- d[, .(catkin_bud_len = mean(catkin_mean, na.rm=T), walnut_diam = mean(walnut_mean, na.rm=T)), by = Variety]

mns <- mns[!is.na(catkin_bud_len)]


# correct names for merge
p[name == 'PI 18256 Manregian', name := 'Manregian']
p[taxa == '59-124', name := '59-124']
p[taxa == '59-165', name := '59-165']
p[taxa == '74-259', name := '74-259']
p[name == 'PI 159568', name := 'PI_159568']
p[name == 'Sir Bon', name := 'SirBon']
p[name == 'Sinensis #5', name := 'Sinensis5']
p[taxa == '56-224', name := '56-224']
p[name == 'Chase D9', name := 'ChaseD9']


mns <- merge(mns, p[, .(Variety = name, protog)])

mns <- rbind(mns[Variety == 'Twister' ,.SD[1]], mns[Variety != 'Twister'])
mns <- rbind(mns[Variety == 'Lara' ,.SD[1]], mns[Variety != 'Lara'])

mns[, protog := ifelse(protog < 0.1, 0, 1)]
mns[, protog := as.factor(protog)]



ggplot(mns, aes(x = protog, y = catkin_bud_len, color = protog)) + 
  geom_jitter(width = 0.1, size = 2) + 
  scale_color_manual(values = c("black", "gray")) + 
  scale_x_discrete(labels = c("protandrous", "protogynous")) + 
  labs(x = '', 
       y = 'Catkin bud length (mm) Jul 20', 
       color = '') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
    legend.position = 'none',
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.y = element_text(size = 14, margin = margin(r = 10))
  )


t.test(x = mns[protog==1, catkin_bud_len],
       y = mns[protog==0, catkin_bud_len])

