library(ggplot2)
library(ggmap)
library(data.table)

library(sp)
library(leaflet)

putah <- fread('~/workspace/heterodichogamy/data/PutahCreek_Jhindsii.tsv')
putah <- putah[, .(ID, type, lat, long)]

putah <- putah[type %in% c("protandrous", "protogynous") & !is.na(long)]

putah1 <- fread("~/Downloads/PutahCreek_Jhindsii - Sheet1.tsv")
putah1 <- putah1[, .(ID, type, lat, long, Prioritize)]
putah1 <- putah1[type %in% c("protandrous", "protogynous") & !is.na(long) & Prioritize == "Y"]


coordinates(putah) <- ~long+lat
leaflet(putah) %>% addMarkers() %>% addTiles()

coordinates(putah1) <- ~long+lat
leaflet(putah1) %>% addMarkers() %>% addTiles()
