library(ggplot2)
library(lme4)
library(data.table)
library(visreg)
library(viridis)

# -----read heterodichogamy phenotypes -----
p <- fread("~/workspace/heterodichogamy/data/heterodichogamy_phenotypes.txt")

# ----- read catkin bud size -----
d <- fread("~/workspace/heterodichogamy/data/New_Stuke_catkin_buds_20Jul2023.tsv")

# ----- read leaf dates -----
l <- fread("~/workspace/heterodichogamy/data/NewStuke_LeafDate_2022.txt")

# ----- format -----
cat_cols <- paste0("catkin", 1:16)
wal_cols <- paste0("walnut", 1:6)

d[, catkin_mean := rowMeans(.SD, na.rm=T), .SDcols = cat_cols]
d[, walnut_mean := rowMeans(.SD, na.rm=T), .SDcols = wal_cols]

mns <- d[, .(catkin_bud_len = mean(catkin_mean, na.rm=T), walnut_diam = mean(walnut_mean, na.rm=T)), by = Variety]

mns <- mns[!is.na(catkin_bud_len)]


# correct names for merge with heterodichogamy phenotype data shet
p[name == 'PI 18256 Manregian', name := 'Manregian']
p[taxa == '59-124', name := '59-124']
p[taxa == '59-165', name := '59-165']
p[taxa == '74-259', name := '74-259']
p[name == 'PI 159568', name := 'PI_159568']
p[name == 'Sir Bon', name := 'SirBon']
p[name == 'Sinensis #5', name := 'Sinensis5']
p[taxa == '56-224', name := '56-224']
p[name == 'Chase D9', name := 'ChaseD9']


# merge with heterodichogamy phenotypes
mns <- merge(mns, p[, .(Variety = name, protog)])

mns <- rbind(mns[Variety == 'Twister' ,.SD[1]], mns[Variety != 'Twister'])
mns <- rbind(mns[Variety == 'Lara' ,.SD[1]], mns[Variety != 'Lara'])

mns[, protog := ifelse(protog < 0.1, 0, 1)]
mns[, protog := as.factor(protog)]


# which names in leaf data sheet need to be corrected for merge?
mns[!Variety %in% l$CULT, Variety]
l[, unique(CULT)]

# correct names for merge
l[CULT == 'Sir Bon', CULT := 'SirBon']
l[CULT == 'Sir Bon', CULT := 'SirBon']
l[CULT == 'PI 18256 Manregian', CULT := 'Manregian']
l[CULT == 'Sinensis #5', CULT := 'Sinensis5']
l[CULT == 'Chase D9', CULT := 'ChaseD9']
l[CULT == 'PI 159568', CULT := 'PI_159568']


mns <- merge(mns, l[, .(Variety = CULT, LFDA)])
mns$LFDA
mns[, lf_date := as.Date(mns$LFDA, format= '%m/%d/%Y')]
mns[, lf := as.numeric(lf_date)]

reference_date <- as.Date("20220301", format = "%Y%m%d")
mns[, days_after_mar1 := as.numeric(difftime(lf_date, reference_date))]


z <- lm(catkin_bud_len ~ days_after_mar1 + protog, mns)
anova(z)
summary(z)
visreg(z)

reference_date <- as.Date("20220301", format = "%Y%m%d")
mns[, days_after_mar1 := as.numeric(difftime(lf_date, reference_date))]

mns[, hist(lf, breaks = 20)]
ggplot(mns, aes(x = protog, y = catkin_bud_len, color = days_after_mar1)) + 
  geom_jitter(width = 0.1, size = 5) + 
  #scale_color_manual(values = c("black", "gray")) + 
  scale_x_discrete(labels = c("protandrous", "protogynous")) + 
  scale_color_viridis() + 
  labs(x = '', 
       y = 'Catkin bud length (mm) Jul 20 2023', 
       color = 'Leafing date\n(days after Mar 1 2022)') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
    #legend.position = 'none',

    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.y = element_text(size = 14, margin = margin(r = 10))
  )


t.test(x = mns[protog==1, catkin_bud_len],
       y = mns[protog==0, catkin_bud_len])





ggplot(mns, aes(x = protog, y = catkin_bud_len)) + 
  geom_jitter(width = 0.1, size = 5) + 
  #scale_color_manual(values = c("black", "gray")) + 
  scale_x_discrete(labels = c("protandrous", "protogynous")) + 
  #scale_color_viridis() + 
  labs(x = '', 
       y = 'Catkin bud length (mm) Jul 20 2023', 
       color = 'Leafing date\n(days after Mar 1 2022)') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        #legend.position = 'none',
        
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14, margin = margin(r = 10))
  )
