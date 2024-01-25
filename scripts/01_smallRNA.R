library(data.table)

sra <- fread("~/workspace/heterodichogamy/smallRNAs/sra.txt", header = F, col.names = c("ID", "tissue"))

cvg <- rbindlist(
  lapply(list.files("~/workspace/heterodichogamy/smallRNAs/results/coverage_JmanNFU/", full.names = T), function(x){
    z <- fread(x, col.names = c("chr", "pos", "cvg"))
    z[, ID := gsub(".txt.gz", "", basename(x))]
  }
  )
)

cvg <- merge(cvg, sra, by = "ID")

genes_Jman_H <- data.table(start = c(32669662,32692488), end = c(32671014,32695292))
genes_Jman_h <- data.table(start = c(31760070,31776222), end = c(31761382,31779026))



# ------ Coverage in windows ------
w <- 1
st <- cvg[, min(pos)]
en <- cvg[, max(pos)]

cvg[, window := cut(pos, breaks = seq(st, en+w, by = w), labels = seq(st, en, by = w), include.lowest =T), by = ID]

cvg_win <- cvg[, .(cvg = mean(cvg)), by = .(ID,window,tissue)]

cvg_win[, window := as.numeric(as.character((window)))]


ggplot(cvg_win, aes(x = window, y = cvg, group = tissue, color = tissue)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(32.67e6, 32.69e6, length.out = 3), labels = sprintf("%.2f",seq(32.67, 32.69, length.out = 3))) +
  labs(color = '', x = 'Position (Mb)', y = "Read depth") +
  scale_color_manual(values = c("#cfb582", "#413a6e")) + 
  theme_classic() + 
  theme(aspect.ratio = .2,
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = c(0.3, 0.9)) + 
  annotate("rect", xmin = genes_Jman_H$start[1],
           xmax = genes_Jman_H$end[1], 
           ymin=80, 
           ymax = 90, 
           alpha = .3,fill = '#94435b') + 
  annotate("rect", xmin = genes_Jman_H$start[2],
           xmax = genes_Jman_H$end[2], 
           ymin=80, 
           ymax = 90, 
           alpha = .3,fill = 'blue') +
  annotate("rect", xmin = 32.6836e6,
           xmax = 32.6916e6, 
           ymin=80, 
           ymax = 90, 
           alpha = .3,fill = 'gray') 
  