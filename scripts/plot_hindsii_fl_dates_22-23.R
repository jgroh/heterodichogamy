library(data.table)
library(readxl)

h <- setDT(read_excel("~/workspace/heterodichogamy/data/PutahCreek_Jhindsii.xlsx"))

h22 <- melt(h, id.vars = c("ID", "type"), measure.vars = c("m_11Apr2022", "f_11Apr2022", "m_18Apr2022", "f_18Apr2022", "m_25Apr2022", "f_25Apr2022"), variable.name = 'date', value.name = 'status')              

h22[, c("sex", "date") := tstrsplit(date, "_")]

h22[, date := as.Date(date, format="%d%b%Y")]

h22 <- h22[type %in% c("protandrous", "protogynous")]

p22 <- h22[, .(proportion = length(which(status == 1))/length(which(is.na(status)))), by = .(date, type, sex)]

p22$interaction <- interaction(p22$type, p22$sex)

p22[interaction == 'protandrous.m', fl := "protandrous_male"]
p22[interaction == 'protandrous.f', fl := "protandrous_female"]
p22[interaction == 'protogynous.m', fl := "protogynous_male"]
p22[interaction == 'protogynous.f', fl := "protogynous_female"]
p22[, fl := factor(fl, levels = c('protandrous_male', 'protandrous_female', 'protogynous_male', 'protogynous_female'))]

ggplot(p22, aes(x = date, y = proportion, fill = fl)) + 
  geom_bar(stat = 'identity') + facet_wrap(~fl, ncol=1)

# -----


h23 <- melt(h, id.vars = c("ID", "type"), measure.vars = c("m_14Apr2023", "f_14Apr2023", "m_21Apr2023", "f_21Apr2023"), variable.name = 'date', value.name = 'status')              

h23[, c("sex", "date") := tstrsplit(date, "_")]

h23[, date := as.Date(date, format="%d%b%Y")]

h23 <- h23[type %in% c("protandrous", "protogynous")]

p23 <- h23[, .(proportion = length(which(status == 1))/length(which(is.na(status)))), by = .(date, type, sex)]

p23$interaction <- interaction(p23$type, p23$sex)

p23[interaction == 'protandrous.m', fl := "protandrous_male"]
p23[interaction == 'protandrous.f', fl := "protandrous_female"]
p23[interaction == 'protogynous.m', fl := "protogynous_male"]
p23[interaction == 'protogynous.f', fl := "protogynous_female"]
p23[, fl := factor(fl, levels = c('protandrous_male', 'protandrous_female', 'protogynous_male', 'protogynous_female'))]

ggplot(p23, aes(x = date, y = proportion, fill = fl)) + 
  geom_bar(stat = 'identity') + facet_wrap(~fl, ncol=1)
