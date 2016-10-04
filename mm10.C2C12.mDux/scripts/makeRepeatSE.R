library(BiocParallel)
library(ShortRead)
library(GenomicAlignments)
BiocParallel::MulticoreParam()
BiocParallel::registered()
library(Rsamtools)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.mDux/inst/scripts/defineParams.R")

pkgDir
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
bamFiles <- bamFiles[grep("C2C12_mDux", bamFiles)]
bamFiles 
bamFiles <- BamFileList(bamFiles)

#'
#' (1) get unique reads with count columns
#' 
library(GenomicAlignments)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/sanitizeReads.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/summarizeOverlaps.adjNH.R")
source("~/tapscott/R_package/repStats/R/repStats.R")
all.reads <- sanitizeReads(bamFiles, worker=6)
save(all.reads, file=file.path(pkgDir, "data", "all.reads.rda"))

#'     
#' (2) count hits: intersectionStrict and adjust for NH
#' 
load("/fh/fast/tapscott_s/CompBio/mm10/mm10.rmsk.rda") ## mm10.rmsk
mm10.rmsk <- keepSeqlevels(mm10.rmsk,
                           value=paste0("chr", c(1:19, "X", "Y")))
features <- getRepNameFeatures(mm10.rmsk)

MinM.repNameSE <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="within",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)

colData(MinM.repNameSE) <-
    DataFrame(BamFile.Path=path(bamFiles),
              NickNames="MinM",
              Treatment=factor(rep(c("DOX", "NODOX"), each=3),
                  levels=c("NODOX", "DOX")))
mcols(MinM.repNameSE) <- as(t(sapply(rowRanges(MinM.repNameSE), .getClassFamily)),
                            "DataFrame")
save(MinM.repNameSE, file=file.path(pkgDir, "data", "MinM.repNameSE.rda"))

## use union mode: type="any"
MinM.repNameUnion <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="any",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)
mcols(MinM.repNameUnion) <- as(t(sapply(rowRanges(MinM.repNameUnion),
                                          .getClassFamily)),
                                 "DataFrame")
colData(MinM.repNameUnion) <- colData(MinM.repNameSE)
save(MinM.repNameUnion, file=file.path(pkgDir, "data", "MinM.repNameUnion.rda"))

#'
#' (3) DESeq and Enrichment Analysis
#'
library(DESeq2)
library(xlsx)
se <- MinM.repNameSE[rowSums(assay(MinM.repNameSE)) > 6, ]
formula <- as.formula("~ Treatment")
ddsRepName <- repDESeq2(se, formula)
sigRepName <- repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
repstats <- repStats(ddsRepName, sigRepName)
a <- summarizeRepStats(repstats)

file <- file.path(pkgDir, "inst", "stats", "MinM.RepeatAnalysis.Beta2.xlsx")
repStatsReport(ddsRepName, sigRepName, repstats,
               file=file, pval=0.05) 

## use counts in union modeR
se <- MinM.repNameUnion[rowSums(assay(MinM.repNameUnion)) > 6, ]
file <- file.path(pkgDir, "inst", "stats", "MinM.RepeatAnalysis.Beta.Union.xlsx")
res <- do.repStats(se, file=file)

#'
#' old RFEA analysis
#' 
library(xlsx)
file <- file.path(pkgDir, "inst", "stats", "MinM.RepeatAnalysis.xlsx")
se <- MinM.repNameSE[rowSums(assay(MinM.repNameSE)) > 6, ]
formula <- as.formula("~ Treatment")
MinM.rfea <- RFEA(se, formula=formula,
                  RFEA.thres=0.1, padj.thres=0.05, verbose=FALSE,
                  file=file)
