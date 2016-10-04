library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")
mm <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")

load("homology.rda")
sum(homology$Gene.ID %in% names(mm)) # 19777
sum(homology$Gene.ID %in% names(hg)) # 18696

humogenes <- c(names(hg), names(mm))
homology <- subset(homology, Gene.ID %in% humogenes)

## do we have one HID to muitiple EntrezID?
tb <- tapply(homology$Gene.ID, homology$HID, function(x) {
             ihg     <- x%in% names(hg)
             imm     <- x %in% names(mm)
             num.hg  <- sum(ihg)
             num.mm  <- sum(imm)
             hgID    <- paste(x[ihg], collapse=";")
             mmID    <- paste(x[imm], collapse=";")
             data.frame(num.hg=num.hg, num.mm=num.mm,hgID=hgID, mmID=mmID)
})

tb <- do.call(rbind, tb)
tb[, "hg:mm"] <- paste(tb$num.hg, tb$num.mm, sep=":")
tb$hgID <- as.character(tb$hgID)
tb$mmID <-  as.character(tb$mmID)

homology.mapping <- tb
save(homology.mapping, file="homology.mapping.rda")
