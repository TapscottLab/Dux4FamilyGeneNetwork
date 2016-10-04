#' must run under R/3.4
library(BiocParallel)
library(ShortRead)
BiocParallel::MulticoreParam(worker=6)
BiocParallel::registered()
library(Rsamtools)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.MB135.HDUX4"
ngsDir <- "/shared/ngs/illumina/jwhiddon/151021_SN367_0569_AH7YMNBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern=".bam$", include.dirs=TRUE,
                       full.names=TRUE,
                       all.files=FALSE)
bamFiles <- bamFiles[1:6]
bamFiles <- BamFileList(bamFiles)

#'
#' Test begines:
#' (1) sanitize the reads: all.reads
#' (2) summarizeOverlaps.2 -> SummarizeExperiments for repeat Names
#' (3) DESeq2
#' 

#' (1) get unique reads with count columns
library(GenomicAlignments)
source("~/tapscott/R_package/repStats/R/repStats.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/sanitizeReads.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/summarizeOverlaps.adjNH.R")
all.reads <- sanitizeReads(bamFiles, worker=6)
save(all.reads, file=file.path(pkgDir, "data", "all.reads.rda"))

#'     
#' (2) count hits: intersectionStrick and adjust for NH
#' 

load("/fh/fast/tapscott_s/CompBio/hg38/hg38.rmsk.rda") ## hg38.rmsk
hg38.rmsk <- keepSeqlevels(hg38.rmsk,
                           value=paste0("chr", c(1:22, "X", "Y")))
##hg38.rmsk <- assignElementID(hg38.rmsk, cores=6)
features <- getRepNameFeatures(hg38.rmsk)

## not adjsut for NH
#tmp <- summarizeOverlaps(features=features,
#                         read=bamFiles,
#                         mode="IntersectionStrict",
#                         inter.feature=TRUE,
#                         ignore.strand=TRUE)

#'
#' colData
#' 
colData <- DataFrame(BamFile.Path=path(bamFiles),
                     NickNames="HinH",
                     Treatment=factor(rep(c("NODOX", "DOX"), 3),
                     levels=c("NODOX", "DOX")))

#'
#' Count hits
#' 
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/summarizeOverlaps.adjNH.R")
## rep element
#HinH.repElementSE <- summarizeOverlaps.adjNH(features=hg38.rmsk,
#                                          all.reads=all.reads,
#                                          type="within",
#                                          inter.feature=TRUE,
#                                             ignore.strand=TRUE)
#colData(HinH.repElementSE) <- colData
#save(HinH.repElementSE, file=file.path(pkgDir, "data", "HinH.repElementSE.rda"))

## repName
HinH.repNameSE <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="within",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)

colData(HinH.repNameSE) <- colData
## a stupid way to get repFamily and repClass
mcols(HinH.repNameSE) <- as(t(sapply(rowRanges(HinH.repNameSE), .getClassFamily)),
                            "DataFrame")
save(HinH.repNameSE, file=file.path(pkgDir, "data", "HinH.repNameSE.rda"))

HinH.repNameUnion <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="any",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)
colData(HinH.repNameUnion) <- colData
mcols(HinH.repNameUnion) <-  as(t(sapply(rowRanges(HinH.repNameUnion),
                                         .getClassFamily)), "DataFrame")
save(HinH.repNameUnion, file=file.path(pkgDir, "data", "HinH.repNameUnion.rda"))

#'
#' (3)DESeq and Enrichment Analysis
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.MB135.HDUX4"
load(file.path(pkgDir, "data", "HinH.repNameSE.rda"))
source("~/tapscott/R_package/repStats/R/repStats.R")
library(DESeq2)
library(xlsx)
se <- HinH.repNameSE[rowSums(assay(HinH.repNameSE)) > 6, ]

formula <- as.formula("~ Treatment")
ddsRepName <- repDESeq2(se, formula)
sigRepName <- repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
repstats <- repStats(ddsRepName, sigRepName)
summarizeRepStats(repstats)
file <- file.path(pkgDir, "inst", "stats", "HinH.RepeatAnalysis.Beta2.xlsx")
repStatsReport(ddsRepName, sigRepName, repstats,
               file=file, pval=0.05)


se <- HinH.repNameUnion[rowSums(assay(HinH.repNameUnion)) > 6, ]
file <- file.path(pkgDir, "inst", "stats", "HinH.RepeatAnalysis.Beta.Union.xlsx")
res <- do.repStats(se, file=file)
                   
#'
#' old codes
#' 
##source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/RFEA.R")
HinH.rfea <- RFEA(se, formula=formula,
                  RFEA.thres=0.1, padj.thres=0.05, verbose=FALSE,
                  file=file)
file <- file.path(pkgDir, "inst", "stats", "HinH.RepeatAnalysis.xlsx")
