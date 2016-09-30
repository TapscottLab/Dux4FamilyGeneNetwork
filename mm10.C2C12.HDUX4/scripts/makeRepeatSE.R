library(BiocParallel)
library(ShortRead)
BiocParallel::MulticoreParam()
BiocParallel::registered()
library(Rsamtools)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4/inst/scripts/defineParams.R")

pkgDir
bamFiles 
bamFiles <- BamFileList(bamFiles)

#' (1) get unique reads with count columns
library(GenomicAlignments)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/sanitizeReads.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/summarizeOverlaps.adjNH.R")
source("~/tapscott/R_package/repStats/R/repStats.R")
all.reads <- sanitizeReads(bamFiles, worker=6)
save(all.reads, file=file.path(pkgDir, "data", "all.reads.rda"))

#'     
#' (2) count hits: intersectionStrick and adjust for NH
#' 
load("/fh/fast/tapscott_s/CompBio/mm10/mm10.rmsk.rda") ## mm10.rmsk
mm10.rmsk <- keepSeqlevels(mm10.rmsk,
                           value=paste0("chr", c(1:19, "X", "Y")))
features <-  getRepNameFeatures(mm10.rmsk)
HinM.repNameSE <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="within",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)

colData(HinM.repNameSE) <-
    DataFrame(BamFile.Path=path(bamFiles),
              NickNames="HinM",
              Treatment=factor(rep(c("DOX", "NODOX"), each=3),
                  levels=c("NODOX", "DOX")))
mcols(HinM.repNameSE) <- as(t(sapply(rowRanges(HinM.repNameSE), .getClassFamily)),
                            "DataFrame")
save(HinM.repNameSE, file=file.path(pkgDir, "data", "HinM.repNameSE.rda"))

## using union mode: type="any"
HinM.repNameUnion <- summarizeOverlaps.adjNH(features=features,
                                          all.reads=all.reads,
                                          type="any",
                                          inter.feature=TRUE,
                                          ignore.strand=TRUE)
colData(HinM.repNameUnion) <- colData(HinM.repNameSE)
mcols(HinM.repNameUnion) <- as(t(sapply(rowRanges(HinM.repNameUnion),
                                        .getClassFamily)),
                            "DataFrame")
save(HinM.repNameUnion, file=file.path(pkgDir, "data", "HinM.repNameUnion.rda"))

#'
#' (3) DESeq and Enrichment Analysis
#'
library(DESeq2)
library(xlsx)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4/inst/scripts/defineParams.R")
source("~/tapscott/R_package/repStats/R/repStats.R")
load(file.path(pkgDir, "data", "HinM.repNameSE.rda"))
load(file.path(pkgDir, "data", "HinM.repNameUnion.rda"))

se <- HinM.repNameSE[rowSums(assay(HinM.repNameSE)) > 6, ]
formula <- as.formula("~ Treatment")
ddsRepName <- repDESeq2(se, formula)
sigRepName <- repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
repstats <- repStats(ddsRepName, sigRepName)
a <- summarizeRepStats(repstats)

file <- file.path(pkgDir, "inst", "stats", "HinM.RepeatAnalysis.Beta.Intersect.xlsx")
repStatsReport(ddsRepName, sigRepName, repstats,
               file=file, pval=0.05) 

## use counts in union modeR
se <- HinM.repNameUnion[rowSums(assay(HinM.repNameUnion)) > 6, ]
file <- file.path(pkgDir, "inst", "stats", "HinM.RepeatAnalysis.Beta.Union.xlsx")
res <- do.repStats(se, file=file)


#'
#' testing
#'
r1 <- read.xlsx(file.path(pkgDir, "inst", "stats", "HinM.RepeatAnalysis.xlsx"),
                sheetIndex=1)
r2 <- read.xlsx(file.path(pkgDir, "inst", "stats", "HinM.RepeatAnalysis.Beta2.xlsx"),
                sheetIndex=1)

#'
#' old codes
#' 
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/RFEA.R")
library(DESeq2)
library(xlsx)
file <- file.path(pkgDir, "inst", "stats", "HinM.RepeatAnalysis.xlsx")

se <- HinM.repNameSE[rowSums(assay(HinM.repNameSE)) > 6, ]
formula <- as.formula("~ Treatment")
HinM.rfea <- RFEA(se, formula=formula,
                  RFEA.thres=0.1, padj.thres=0.05, verbose=FALSE,
                  file=file)
