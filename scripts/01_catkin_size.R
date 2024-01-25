library(data.table)
library(ggplot2)
jh <- fread("~/workspace/heterodichogamy/data/PutahCreek_Jhindsii.tsv")


# first look
ggplot(data = jh[type %in% c("protandrous","protogynous")], 
       aes(x = L1_21Apr2023, 
           y = W1_21Apr2023, 
           color = type)) + geom_point()

#jh[, catkin_len := rowMeans(.SD[, .(L1_14Apr2023, L2_14Apr2023, L3_14Apr2023)], na.rm=T)]
#jh[, catkin_wd := rowMeans(.SD[, .(W1_14Apr2023, W2_14Apr2023, W3_14Apr2023)], na.rm=T)]


ggplot(jh[type %in% c("protandrous", "protogynous")], 
       aes(x = L1_21Apr2023, y = W1_21Apr2023, shape = type, color = type)) +
         geom_point(size = 2) + 
  scale_color_manual(values = c("#cfb582", "#413a6e")) +
  theme_classic() + 
  labs(x = 'Catkin length Apr 21 (mm)',
       y = 'Catkin width Apr 21 (mm)', 
       shape = '') + 
  #scale_y_continuous(breaks = c(0,10,15,20)) +
  theme(aspect.ratio = 1, 
        plot.margin = margin(30,30,30,30, "pt"),
        text = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2),
        #legend.position = c(0.15, 0.85), 
        legend.position = 'none')

# pca
pc_data <- na.omit(jh[, .(type,catkin_len, catkin_wd)])
PCvals <- prcomp(pc_data[, .(catkin_len, catkin_wd)])$x
PC <- data.table(type = pc_data[, type], PCvals)
ggplot(PC[type %in% c("protandrous", "protogynous")], 
       aes(x = PC1, y = PC2, color = type)) + geom_point()

lda
