## This script read the bam files and create a SummarizedExperiments instance containing
## the genomic position and count hits.

library(GenomicFeatures)
library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicAlignments)
library(BiocParallel)
library(ShortRead)
BiocParallel::MulticoreParam()
BiocParallel::registered()

source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.mDux/inst/scripts/defineParams.R")

bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
bamFiles <- bamFiles[grep("C2C12_mDux", bamFiles)]

###
## Make SummarizedExprients with count reads
###    
genes.mm10 <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
MinM.knownGene <- summarizeOverlaps(features=genes.mm10,
                                      reads=BamFileList(bamFiles),
                                      mode="IntersectionStrict",
                                      inter.feature=TRUE,
                                      singleEnd=TRUE,
                                      ignore.strand=TRUE) 


sampleInfo <- data.frame(sample_name = sub(".bam", "", names(bamFiles)),
                     file_bam = path(bamFiles),
                     nick_name = "MinM",
                     Treatment=factor(rep(c("DOX", "NODOX"), each=3),
                          levels=c("NODOX", "DOX")),
                     stringsAsFactors=FALSE)

si <- SGSeq::getBamInfo(sampleInfo, core=6)
save(si, file=file.path(pkgDir, "data", "si.rda"))
colData(MinM.knownGene) <- as(si, "DataFrame")

source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/getScaledCounts.R")
assays(MinM.knownGene)$FPKM <- getScaledCountsPerTx(MinM.knownGene)

save(MinM.knownGene, file=file.path(pkgDir, "data", "MinM.knownGene.rda"))

## print out FPKM for Jenn, along with Symbol
id <- names(rowRanges(MinM.knownGene))
sym <- select(org.Mm.eg.db, keys=id, keytype="ENTREZID", columns="SYMBOL")
output <- cbind(sym, as.data.frame(assays(MinM.knownGene)$FPKM, "data.frame"))
write.csv(output, file=file.path(pkgDir, "inst", "extdata", "SE.FPKM.csv"))
          
