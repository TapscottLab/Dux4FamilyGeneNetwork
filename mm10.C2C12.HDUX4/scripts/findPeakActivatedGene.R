#!/usr/bin/env Rscript
#SBATCH -N1 -n12 -t 1-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s
#./findPeakActivatedGene.R

source("defineParams.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/findAltPromoter.R")
library(Rsamtools)
library(GenomicAlignments)
cores <- 12
load("/fh/fast/tapscott_s/CompBio/ChIP-Seq/mm10.HuMoDoxy2/data/peaks_grl.rda")
n <- "Sample_2-HinM-WITHdoxy-MO-24:Sample_1-HinM-WITHdoxy-HA-24_peaks"
peaks <- peaks_grl[[n]]
names(peaks) <- paste0("peak_", 1:length(peaks))

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
exons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
tx <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")

#'
#' define bins 
#' 
si <- seqinfo(Rsamtools::BamFileList(bamFiles))
sl <- seqlevels(si)
which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), strand="*")
bins <- split(which, seq_along(which))

#'
#' to linkPeaksToTx: get reads links peaks to transcripts
#'
links <- lapply(bamFiles, linkPeaksToTx, list_which=bins,
                peaks=peaks, tx=exons, mc.cores=1)
names(links) <- sub(".bam", "", basename(bamFiles))
splicing_links <- lapply(bamFiles, splicePeaksToTx, list_which=bins,
                         peaks=peaks, tx=tx, mc.cores=1)
names(splicing_links) <- sub(".bam", "", basename(bamFiles))

#'
#' get hits of the common links; return hits of the peak-gene pair
#' 
link_hits <- getCommonLinksHits(dox_links=links[1:3], nodox_links=links[4:6])

splicing_hits <- getCommonLinksHits(dox_links=splicing_links[1:3],
                                    nodox_links=links[4:6])

#'
#' evaluate activation score and annotatoin
#'
library(org.Mm.eg.db)
df <- linkAnnotation(org=org.Mm.eg.db, link_hits, splicing_hits, pkgDir)

################## end of pipeline ####################


 
