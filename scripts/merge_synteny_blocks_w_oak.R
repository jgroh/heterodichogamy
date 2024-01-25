# ----- 3 way synteny w oak ----- 
AB <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Qlobata.Lakota_v1.i2.blocks.2', 
            header=F, sep = '\t',
            col.names = c("A", "B"))


BC <- fread('~/workspace/heterodichogamy/HC_locus_structure/MCscan/Lakota_v1.Pawnee.i1.blocks', 
            header=F, sep = '\t',
            col.names = c("B", "C"))

ABC.blocks <- merge(AB, BC, all= T, by = 'B')

ABC.blocks[is.na(A), A := '.']
ABC.blocks[is.na(B), B := '.']
ABC.blocks[is.na(C), C := '.']

fwrite(ABC.blocks[, .(A, B, C)], 
       file = '~/workspace/heterodichogamy/HC_locus_structure/MCscan/Oak.Lak1.Pawnee.blocks.all.manual',
       col.names = F,
       row.names = F, 
       quote = F,
       sep = '\t')

