library(data.table)
library(ggplot2)


wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/wgs_samples.tsv")
gbs_samples <- fread("~/workspace/heterodichogamy/pecan/GBS_data/metadata.tsv")
kinship1 <- fread("~/workspace/heterodichogamy/pecan/WGS_GBS_samples.relatedness")
kinship2 <- fread("~/workspace/heterodichogamy/pecan/WGS_GBS_samples.relatedness2")

# inspect intersection of varieties
sort(wgs_samples[, type])
sort(gbs_samples[, Cultivar])
sort(intersect(wgs_samples[, type], gbs_samples[, Cultivar]))

# the matrices have redundant rows. We are expecting 34*83 unique combinations
# Note matrices don't have same format
kinship1 <- kinship1[INDV1 %in% gbs_samples[, Accession] & INDV2 %in% wgs_samples[, run]]

# the 2nd kinship matrix has redundant rows
kinship2 <- kinship2[INDV1 %in% wgs_samples[, run] & INDV2 %in% gbs_samples[, Accession]]

# Assign variety in 1st kinship matrix 
for(i in 1:nrow(kinship1)){
  ind1 <- kinship1[i, INDV1]
  ind2 <- kinship1[i, INDV2]
  
  kinship1[i, variety1 := gbs_samples[Accession == ind1, Cultivar]]
  kinship1[i, variety2 := wgs_samples[run == ind2, type]][]
}
kinship1

# Assign variety in 2nd kinship matrix 
for(i in 1:nrow(kinship2)){
  ind1 <- kinship2[i, INDV1]
  ind2 <- kinship2[i, INDV2]
  
  kinship2[i, variety1 := wgs_samples[run == ind1, type]]
  kinship2[i, variety2 := gbs_samples[Accession == ind2, Cultivar]]
}


# assign variable to indicate whether labelled as same variety
kinship1[variety1 != variety2, sm := 0]
kinship1[variety1 == variety2, sm := 1]

kinship2[variety1 != variety2, sm := 0]
kinship2[variety1 == variety2, sm := 1]


# look at relatedness statistics for samples labelled as same variety vs. not
ggplot(kinship1[variety2 == 'PoSey'], aes(y = RELATEDNESS_AJK, x = as.factor(sm))) + 
  geom_point(position = position_jitter(width = 0.2))

# kinship1[variety2 == 'PoSey'] The other homozygote appears to just be a clone of Mahan

# these individuals look ok to use
keep1 <- kinship1[sm == 1 & RELATEDNESS_AJK > 0.6, variety1]


kinship1[sm == 0 & RELATEDNESS_AJK > 0.6]

kinship2[, k1 := ifelse(variety1 %in% keep1 & variety1 == variety2, 'red', 'black')]

kinship2[, k1 := ifelse(variety1 %in% keep1 & variety1 == variety2, 'red', 'black')]

ggplot(kinship2, aes(y = RELATEDNESS_PHI, x = as.factor(sm), color = k1)) + 
  geom_point(position = position_jitter(width = 0.2))






