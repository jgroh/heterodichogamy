library(data.table)

library(data.table)
library(ggplot2)
library(cowplot)


Jcal_pstart <- 31348000
Jcal_pend <- 31400000


# ===== Read coverage =====


cvg <- rbindlist(lapply(list.files("~/workspace/heterodichogamy/smallRNAs/results/coverage_Jcali_primary/", 
                                                  pattern = '*.txt.gz', full.names = T),
                                       function(x){
                                         z <- fread(x, select = 2:3, col.names = c("position", "coverage"))
                                         z[, sample := gsub(".txt.gz", "", basename(x))]
                                         return(z)
                                       }))


# ---- calculate coverage in 10bp windows -----
cvg[, window10 := cut(position, breaks = seq(Jcal_pstart, Jcal_pend+10, by = 10), labels = seq(Jcal_pstart, Jcal_pend, by = 10), include.lowest =T), by = sample]

cvg10bp <- cvg[, .(coverage = mean(coverage)), by = .(sample, window10)]

cvg10bp[, window10 := as.numeric(as.character((window10)))]


ggplot(cvg10bp[coverage != 0], aes(x = window10, y = coverage, group = sample, color = sample)) + geom_point() 

ggplot(cvg10bp[coverage != 0 & sample == 'walnut_smRNA_JP3'], aes(x = window10, y = coverage, group = sample, color = sample)) + geom_point() 

cvg10bp[coverage != 0 & sample == 'walnut_smRNA_JP3']

ggplot(cvg10bp[coverage != 0 & window10 > 31360000 & window10 < 31380000], aes(x = window10, y = coverage, 
                                                               group = sample, color = sample)) + geom_point() + geom_line()



ggplot(cvg10bp[coverage != 0 & window10 > 31360000 & window10 < 31380000 & sample == 'walnut_smRNA_JP3'], aes(x = window10, y = coverage, 
                                                                               group = sample, color = sample)) + geom_point() + geom_line()


ggplot(cvg[coverage !=0 & position > 31360000 & position < 31380000], aes(x = position, y = coverage, 
                                                               group = sample, color = sample)) + geom_point() + geom_line()


ggplot(cvg[coverage !=0 & position > 31360000 & position < 31365000], aes(x = position, y = coverage, 
                                                                          group = sample, color = sample)) + geom_point() + geom_line()

