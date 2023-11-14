f1 <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Qlobata.Lakota_v1.i2.blocks.manual', 
            header=F, sep = '\t',
            col.names = c("Q", "L"))
f2 <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Pawnee.i1.blocks', 
            header=F, sep = '\t',
            col.names = c("L", "P"))

blocks <- merge(f1, f2, all= T, by = 'L')

blocks[is.na(Q), Q := '.']
blocks[is.na(P), P := '.']

fwrite(blocks[, .(Q, L, P)], 
       file = '~/workspace/heterodichogamy/HC_locus_structure/MCscan/3way.blocks.manual.edit',
       col.names = F,
       row.names = F, 
       quote = F,
       sep = '\t')
