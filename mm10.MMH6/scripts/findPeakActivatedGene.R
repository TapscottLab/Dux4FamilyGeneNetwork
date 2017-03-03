#!/app/easybuild/software/R/3.3.1-foss-2016b-fh1/bin/Rscript
#./findPeakActivatedGene.R
#SBATCH -N1 -n12 -t 2-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s


#' module load R/3.3.1-foss-2016b-fh1

source("~/tapscott/RNA-Seq/R_scripts/findAltPromoter.R")
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(Rsamtools)
library(GenomicAlignments)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
cores <- 12

#' define bamFiles
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
ngsDir <- "/shared/ngs/illumina/jwhiddon/161121_D00300_0344_AH5MTNBCXY"
bamDir <- file.path(ngsDir, "tophat", "bam")
MMH_bamFiles <- list.files(bamDir, pattern="\\.bam$",  include.dirs=TRUE,
                       full.names=TRUE)

ngsDir <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
mDux_bamFiles <- bamFiles[grep("C2C12_mDux", bamFiles)]
luc_bamFiles <- bamFiles[grep("Luc", bamFiles)]

#' define features
tx <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
exons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")

#' define peaks
load("/fh/fast/tapscott_s/CompBio/DownloadSeqData/MMH6/data/peaks_grl.rda")
peaks <- peaks_grl[[1]]
names(peaks) <- paste0("peak_", 1:length(peaks))

#'
#' define bins 
#' 
si <- seqinfo(Rsamtools::BamFileList(MMH_bamFiles))
sl <- seqlevels(si)
which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), strand="*")
bins <- split(which, seq_along(which))

#'
#' find Peak Activated Genes pairs
#' 

## after change the minoverlap=2L, this works fine.
findPeakActivatedGene <- function(bamFiles, bins, peaks, cores, tx, exons) {
    links <- lapply(bamFiles, linkPeaksToTx, list_which=bins,
                    peaks=peaks, tx=exons, mc.cores=cores)
    names(links) <- sub(".bam", "", basename(bamFiles))
    splicing_links <- lapply(bamFiles, splicePeaksToTx, list_which=bins,
                             peaks=peaks, tx=tx, mc.cores=cores)
    names(splicing_links) <- sub(".bam", "", basename(bamFiles))
    link_hits <- getCommonLinksHits(dox_links=links[4:6], nodox_links=links[1:3])
    splicing_hits <- getCommonLinksHits(dox_links=splicing_links[4:6],
                                    nodox_links=links[1:3])
    library(org.Mm.eg.db)
    df <- linkAnnotation(org=org.Mm.eg.db, link_hits, splicing_hits, pkgDir, peaks)
}

## MMH6
bamFiles <- c(luc_bamFiles[4:6], MMH_bamFiles[4:6]) ## first three be nodox and last three are dox
df <- findPeakActivatedGene(bamFiles, bins=bins, peaks=peaks, cores=cores, tx=tx, exons=exons)
cmd <- paste("mv ", file.path(pkgDir, "inst", "stats", "PeakActivatedGene.csv"),
              file.path(pkgDir, "inst", "stats", "MMH6_PeakActivatedGene.csv"))
system(cmd)
cmd <- paste("mv ", file.path(pkgDir, "data", "PeakActivatedGene.rda"),
              file.path(pkgDir, "data", "MMH6_PeakActivatedGene.rda"))
system(cmd)

## MinM
bamFiles <- c(luc_bamFiles[4:6], mDux_bamFiles[1:3]) ## first three be nodox and last three are dox
df <- findPeakActivatedGene(bamFiles, bins, peaks, cores, tx=tx, exons=exons)
cmt <- paste("mv ", file.path(pkgDir, "inst", "stats", "PeakActivatedGene.csv"),
              file.path(pkgDir, "inst", "stats", "MinM_PeakActivatedGene.csv"))
system(cmt)
cmt <- paste("mv ", file.path(pkgDir, "data", "PeakActivatedGene.rda"),
              file.path(pkgDir, "data", "MinM_PeakActivatedGene.rda"))
system(cmt)

#'
#' a test for the bug in pinterest: reproducing the error
#' 
#bam_file="/shared/ngs/illumina/jwhiddon/161121_D00300_0344_AH5MTNBCXY/tophat/bam/6_MMH6_WITHdox_rep3.bam"
#'x=bins[[19]]
#'flag <- scanBamFlag(isSecondaryAlignment = FALSE)
#'param <- ScanBamParam(flag = flag, tag = "XS", which = x)
#'gap <- GenomicAlignments::readGAlignments(file = bam_file,
#'                                          param = param)
#'ol <- .bestOverlaps(query=gap, subject=peaks, ignore.strand=TRUE,
#'                    minoverlap=1L)
#'query=gap
#'pinterest(query[359908], subject[47658])
