## This script read the bam files and create a SummarizedExperiments instance containing
## the genomic position and count hits.

library(GenomicFeatures)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicAlignments)
library(BiocParallel)
library(ShortRead)
BiocParallel::MulticoreParam()
BiocParallel::registered()

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.MB135.HDUX4"
ngsDir <- "/shared/ngs/illumina/jwhiddon/151021_SN367_0569_AH7YMNBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern=".bam$", include.dirs=TRUE,
                       full.names=TRUE,
                       all.files=FALSE)
bamFiles <- bamFiles[1:6]
bamFiles <- BamFileList(bamFiles)

###
## Make SummarizedExprients with count reads
###    
genes.hg38 <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")
HinH.knownGene <- summarizeOverlaps(features=genes.hg38,
                                    reads=bamFiles,
                                    mode="IntersectionStrict",
                                    inter.feature=TRUE,
                                    singleEnd=TRUE,
                                    ignore.strand=TRUE) 

sampleInfo <- data.frame(sample_name = sub(".bam", "", names(bamFiles)),
                         file_bam=path(bamFiles),
                         nick_name="HinH",
                         Treatment=factor(rep(c("NODOX", "DOX"), 3),
                             levels=c("NODOX", "DOX")),
                         stringsAsFactors=FALSE)
library(SGSeq)
si <- getBamInfo(sampleInfo, core=6)
save(si, file=file.path(pkgDir, "data", "si.rda"))
colData(HinH.knownGene) <- as(si, "DataFrame")

source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/getScaledCounts.R")
assays(HinH.knownGene)$FPKM <- getScaledCountsPerTx(HinH.knownGene)

#' validation of FPKM - done. looks great.
save(HinH.knownGene, file=file.path(pkgDir, "data", "HinH.knownGene.rda"))
## print out FPKM for Jenn, along with Symbol
id <- names(rowRanges(HinH.knownGene))
sym <- select(org.Hs.eg.db, keys=id, keytype="ENTREZID", columns="SYMBOL")
output <- cbind(sym, as.data.frame(assays(HinH.knownGene)$FPKM, "data.frame"))
write.csv(output, file=file.path(pkgDir, "inst", "extdata", "SE.FPKM.csv"))

###
## RMSK SE
###
load("/fh/fast/tapscott_s/CompBio/hg38/hg38.rmsk.rda")
SE.rmsk <- summarizeOverlaps(features=hg38.rmsk,
                             reads=bamFiles,
                             mode="IntersectionStrict",
                             inter.feature=TRUE,
                             singleEnd=TRUE,
                             ignore.strand=TRUE)
sampleNames <- sub(".bam", "", colnames(SE.rmsk))
colData <- DataFrame(BamFile.Path=path(bamFiles),
                     Sample.Names=sampleNames,
                     NickNames="HinH",
                     Treatment=factor(rep(c("NODOX", "DOX"), 3),
                         levels=c("NODOX", "DOX")))
colData(SE.rmsk) <- colData
save(SE.rmsk, file=file.path(pkgDir, "data", "SE.rmsk.rda"))

###
## RMSK repName: hit count should be adjusted by NH (write another summarizeOverlaps)
###
features <- split(hg38.rmsk, hg38.rmsk$repName)

SE.repName <- summarizeOverlaps(features=features,
                                reads=bamFiles,
                                mode="IntersectionStrict",
                                inter.feature=TRUE,
                                singleEnd=TRUE,
                                ignore.strand=TRUE)

features.repFam <- split(hg38.rmsk, hg38.rmsk$repFamily)
SE.repFamily <- summarizeOverlaps(features=features.repFam,
                                  read=bamFiles,
                                  mode="IntersectionStrict",
                                  inter.feature=TRUE,
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE)
colData(SE.repName) <- colData(SE.repFamily) <- colData
save(SE.repName, file=file.path(pkgDir, "data", "SE.repName.rda"))
save(SE.repFamily, file=file.path(pkgDir, "data", "SE.repFamily.rda"))



