#!/usr/bin/env Rscript
#SBATCH -N1 -n12 -t 1-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s
#./findPeakActivatedGene.R

source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.mDux/inst/scripts/defineParams.R")
source("~/tapscott/RNA-Seq/R_scripts/findAltPromoter.R")
library(Rsamtools)
library(GenomicAlignments)

cores <- 12
load("/fh/fast/tapscott_s/CompBio/ChIP-Seq/mm10.HuMoDoxy/data/peaks_grl.rda")
n <- "AllSWAP_DOXY_vs_AllControl_peaks"
peaks <- peaks_grl[[n]]
names(peaks) <- paste0("peak_", 1:length(peaks))

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
tx <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
exons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")

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


#' (4) visualization for df
#'source("~/tapscott/RNA-Seq/R_scripts/findAltFirstExon.R")
#'df_LTR <- df[df$Peak.ov.LTR, ]
#'figDir <- file.path(pkgDir, "inst", "figures", "Peak_Activated_Gene_LTR")
#'mclapply(rownames(df_LTR), function(x) {
#'    symbol=df_LTR[x, "Gene.Symbol"]
#'    pdf(file.path(figDir, paste0(symbol, ":", x, ".pdf")))
#'    plotFeatureBySymbol(symbol, sgfc, org=org.Mm.eg.db)
#'    dev.off()
#'}, mc.cores=12, mc.preschedule=FALSE)
        
#'sub <- df[df$Peak.ov.LTR & df$Splicing.link, ]
#'pdf(file.path(figDir, "splicing_link.pdf"))
#'lapply(rownames(sub), function(x) {
#'    symbol=df_LTR[x, "Gene.Symbol"]
#'    plotFeatureBySymbol(symbol, sgfc, org=org.Mm.eg.db)
#'})
#'dev.off()


#'
#' compare to SGSeq
#'
#'sgseq <- get(load(file.path(pkgDir, "data", "sgvc.AS.LTR.report.rda")))
#'sgseq <- sgseq$E_df
#'df_LTR <- df[df$Peak.ov.LTR, ]
#'sum(df_LTR$Gene.EntrezID %in% as.character(sgseq$EntrezID))

#'sdf_LTR <- splicing_df[splicing_df$Peak.ov.LTR, ]
#'sum(splicing_df$Gene.EntrezID %in% as.character(sgseq$EntrezID))
              
