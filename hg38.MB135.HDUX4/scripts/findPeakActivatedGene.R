#!/usr/bin/env Rscript
#SBATCH -N1 -n12 -t 1-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s
#./findPeakActivatedGene.R
library(Rsamtools)
library(GenomicAlignments)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.MB135.HDUX4"
ngsDir <- "/shared/ngs/illumina/jwhiddon/151021_SN367_0569_AH7YMNBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern=".bam$", include.dirs=TRUE,
                       full.names=TRUE,
                       all.files=FALSE)
bamFiles <- bamFiles[1:6]


source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/findAltPromoter.R")
cores <- 12
load("/fh/fast/tapscott_s/CompBio/ChIP-Seq/hg38.MB135.dux4/MACS/hMB1354DFLvshMB135Dneg.rda") ## peaks

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")
tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")

#'
#' define bins 
#' 
si <- seqinfo(Rsamtools::BamFileList(bamFiles))
sl <- seqlevels(si)[1:24]
which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), strand="*")
bins <- split(which, seq_along(which))

#'
#' to linkPeaksToTx: get reads links peaks to transcripts
#'
names(peaks) <- peaks$name
links <- lapply(bamFiles, linkPeaksToTx, list_which=bins,
                peaks=peaks, tx=exons, mc.cores=12)
names(links) <- sub(".bam", "", basename(bamFiles))
splicing_links <- lapply(bamFiles, splicePeaksToTx, list_which=bins,
                         peaks=peaks, tx=tx, mc.cores=12)

names(splicing_links) <- sub(".bam", "", basename(bamFiles))
#'
#' get hits of the common links; return hits of the peak-gene pair
#' 
link_hits <- getCommonLinksHits(dox_links=links[c(2,4,6)], nodox_links=links[c(1,3,5)])

splicing_hits <- getCommonLinksHits(dox_links=splicing_links[c(2,4,6)],
                                    nodox_links=splicing_links[c(1,3,5)])
## splicing hits has 296

#'
#' evaluate activation score and annotatoin
#'
library(org.Hs.eg.db)
df <- linkAnnotation(org=org.Hs.eg.db, link_hits=link_hits,
                     splicing_hits=splicing_hits, pkgDir=pkgDir,
                     peaks=peaks)

df_splicing <- splicingLinkAnnotation(org=org.Hs.eg.db,
                                      link_hits=link_hits,
                                      splicing_hits=splicing_hits,
                                      pkgDir=pkgDir, peaks=peaks)

################## end of pipeline ####################

#'
#' Compare to Zizhen's list
#' 
library(xlsx)
file <- file.path(pkgDir, "inst", "extdata", "journal.pgen.1003947.s017.xlsx")
s4 <- read.xlsx(file, sheetIndex=1, startRow=2)
sum(df_splicing$Gene.Symbol %in% s4$gene.symbol)
sum(df$Gene.Symbol %in% s4$gene.symbol)

i <- which(!df$Gene.Symbol %in% s4$gene.symbol)
j <-  which(!df_splicing$Gene.Symbol %in% s4$gene.symbol)
## checked chr19:7,030,404-7,030,610
## checked hMB1354DFLvshMB135Dneg_peak_7528-2962 (34)
##         peak range - chr19:6,390,818-6,393,328 (need to verify) - GTF2F1
##         UCSC genome - chr19:6,390,818-6,393,328
## need to ignore MIRNA
