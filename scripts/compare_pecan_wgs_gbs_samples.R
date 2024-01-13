library(data.table)
library(ggplot2)

# ---- read data -----
wgs_samples <- fread("~/workspace/heterodichogamy/pecan/WGS_data/all_wgs_varieties_unverified.txt")
gbs_samples <- fread("~/workspace/heterodichogamy/pecan/GBS_data/metadata.tsv")
analysis_set <- fread("~/workspace/heterodichogamy/pecan/WGS_data/analysis_set.txt")
kinship1 <- fread("~/workspace/heterodichogamy/pecan/WGS_GBS_samples.relatedness")
kinship2 <- fread("~/workspace/heterodichogamy/pecan/WGS_GBS_samples.relatedness2")

# ---- clean data -----
# Note matrices don't have same format
# there are 53+83=136 individuals
# first matrix has choose(136,2) + 136 = 9316 rows

# the 2nd kinship matrix has redundant rows
# the 2nd matrix has n^2 rows, where n is number of individuals in vcf
# filter these to remove duplicates

keep_index <- which(kinship2[, paste0(INDV1, INDV2, sep = '_')] %in% kinship1[, paste0(INDV1, INDV2, sep = '_')])
kinship2 <- kinship2[keep_index]

kinship <- merge(kinship1, kinship2)


# ----- assign varieties -----

# Assign varieties in kinship matrix 
for(i in 1:nrow(kinship)){
  ind1 <- kinship[i, INDV1]
  ind2 <- kinship[i, INDV2]
  
  if(ind1 %in% gbs_samples[, Accession]){
    kinship[i, source1 := 'GBS_SRA']
    kinship[i, variety1 := gbs_samples[Accession == ind1, Cultivar]]
  }
  if(ind2 %in% gbs_samples[, Accession]){
    kinship[i, source2 := 'GBS']
    kinship[i, variety2 := gbs_samples[Accession == ind2, Cultivar]]
  }
  
  if(ind1 %in% wgs_samples[, ID]){
    kinship[i, source1 := 'WGS']
    kinship[i, variety1 := wgs_samples[ID == ind1, variety]]
  }
  if(ind2 %in% wgs_samples[, ID]){
    kinship[i, source2 := 'WGS']
    kinship[i, variety2 := wgs_samples[ID == ind2, variety]]
  }
  
}
kinship

# ---- assign identity status ----

# set variable to indicate whether two individuals are same sample in the vcf
kinship[variety1 != variety2, identity := 0]
kinship[variety1 == variety2 & INDV1 != INDV2, identity := 1]
kinship[INDV1==INDV2, identity := 2]



# ---- plots -----
ggplot(kinship, aes(x = as.factor(identity), y = RELATEDNESS_PHI)) + geom_boxplot()
ggplot(kinship, aes(x = as.factor(identity), y = RELATEDNESS_AJK)) + geom_point()


newdata <- kinship[(grepl('CILL', INDV1) & !grepl("CILL", INDV2)) | (!grepl("CILL", INDV1) & grepl("CILL", INDV2)) | identity == 2]
ggplot(newdata, aes(x = as.factor(identity), y = RELATEDNESS_PHI)) + geom_point()
ggplot(newdata, aes(x = as.factor(identity), y = RELATEDNESS_AJK)) + geom_point()


newdata[identity == 1 & RELATEDNESS_AJK < 0.4]
newdata[identity == 1 & RELATEDNESS_AJK > 0.4]

newdata[identity == 0 & RELATEDNESS_PHI > 0.4]
newdata[identity == 0 & RELATEDNESS_AJK > 0.75]
newdata[identity == 0 & RELATEDNESS_AJK > 0.5]


# ----- subset to analysis set -----
kinship_sub <- kinship[INDV1 %in% analysis_set[, ID] & INDV2 %in% analysis_set[, ID]]
ggplot(kinship_sub, aes(x = as.factor(identity), y = RELATEDNESS_AJK)) + geom_point()
ggplot(kinship_sub, aes(x = as.factor(identity), y = RELATEDNESS_PHI)) + geom_jitter() 


kinship_sub[identity == 0 & RELATEDNESS_AJK > 0.6] # will drop Major from new sequencing when doing GWAS


kinship_sub[identity == 0 & RELATEDNESS_PHI > 0.2] # will drop Major from new sequencing when doing GWAS

kinship_sub[variety1 == 'Mahan'] # will drop Major from new sequencing when doing GWAS



# look at relatedness statistics for samples labelled as same variety vs. not
ggplot(kinship1, aes(y = RELATEDNESS_AJK, x = as.factor(sm))) + 
  geom_point(position = position_jitter(width = 0.2))

# kinship1[variety2 == 'PoSey'] The other homozygote appears to just be a clone of Mahan

# these individuals look ok to use
keep1 <- kinship1[sm == 1 & RELATEDNESS_AJK > 0.6, variety1]


kinship1[sm == 0 & RELATEDNESS_AJK > 0.6]

kinship2[, k1 := ifelse(variety1 %in% keep1 & variety1 == variety2, 'red', 'black')]

kinship2[, k1 := ifelse(variety1 %in% keep1 & variety1 == variety2, 'red', 'black')]

ggplot(kinship2, aes(y = RELATEDNESS_PHI, x = as.factor(sm), color = k1)) + 
  geom_point(position = position_jitter(width = 0.2))






