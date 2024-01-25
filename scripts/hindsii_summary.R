library(readxl)
library(data.table)

p <- fread("~/workspace/heterodichogamy/data/PutahCreek_Jhindsii.tsv")

p[type %in% c('protogynous', 'protandrous'), .N, by = type]

h <- setDT(read_excel("~/workspace/heterodichogamy/Manuscript/Supplementary_tables.xlsx", sheet = 'NovaSeq'))


h[species == 'hindsii', .N, by = phenotype]
