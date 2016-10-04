#'
#' This script mimics Zizhen's algorithm in Janet Young's paper
#' NOTE: it is only for single-ended reads
#' 
#' This script aims to find genes activated by peaks. A peak-gene list is yield by
#' using our ChIP-seq and RNA-seq data. The algorithm follows the steps below:
#' linkPeaksToTx():
#' (1) retain reads overlaps with peaks
#' (2) retain those reads overlaps with genes (exons) with at least 20 bps 
#' (3) create peak-gene list using those reads
#'
#' An alternative step is to find the reads splicing into transcripts (indel or exons)
#' and then creat the peak-gene list (splicePeaksToTx)
#'

.removeDupMapping <- function(ol) {
    ## must be 1:1 (query:subject) mapping. If 1:n (n>1), chose the best one
    dup <- duplicated(queryHits(ol))
    ol <- ol[!dup, ]
    ol
}

.bestOverlaps <- function(query, subject, ignore.strand=TRUE, minoverlap=1L) {
    ## We want one query hits ony one subject. If there is 1-to-n mapping
    ## where n > 1, the function picks the best subject based upon the overlapping
    ## ranges. The function yields n-to-1 mapping, where n >= 1.
    ## NOTE: both query and subject have to be GRanges, not List-type.
    ol <- findOverlaps(query, subject, ignore.strand=ignore.strand,
                       minoverlap=minoverlap)
    w <- width(pintersect(query[queryHits(ol)], subject[subjectHits(ol)]))
    best <- tapply(w, queryHits(ol), function(x) {
        keep <- rep(FALSE, length(x))
        keep[which.max(x)[1]] <- TRUE
        keep
    })
    best <- as.logical(unlist(best))
    ol.best <- ol[best]
}

## test: ol[which(queryHits(ol)==2323)]
.bestOverlaps_grl <- function(query, subject, ignore.strand=TRUE, minoverlap=1L) {
    ## We want one query hits ony one subject. If there is 1-to-n mapping
    ## where n > 1, the function picks the best subject based upon the overlapping
    ## ranges. The function yields n-to-1 mapping, where n >= 1.
    ## NOTE: both query and subject have to be GRanges, not List-type.
    ol <- findOverlaps(query, subject, ignore.strand=ignore.strand,
                       minoverlap=minoverlap)
    ## which reads hit multiple genes (exons)? remove multiple hits
    if (any(duplicated(queryHits(ol)))) {
        tb <- table(queryHits(ol))
        mul_hits <- names(tb)[tb > 1]
        not_keep <- sapply(mul_hits, function(mhits) {
            idx <- which(queryHits(ol) == as.integer(mhits))
            x <- query[queryHits(ol)[idx[1]]]
            pin <- lapply(subject[subjectHits(ol)[idx]], pintersect, as(x, "GRanges"))
            w <- sapply(pin, function(p) sum(width(p)))
            idx[which.min(w)[1]]
        })
        ol <- ol[-not_keep]
    }
    ol
    
}

splicePeaksToTx <- function(bam_file, list_which, peaks, tx,
                            mc.cores=1) {
    message(basename(bam_file))
    ## this is for dox samples and find the peak-gene pairs
    list_link_reads <- mclapply(list_which, function(x) {
        flag <- scanBamFlag(isSecondaryAlignment = FALSE)
        param <- ScanBamParam(flag = flag, tag = "XS", which = x)
        gap <- GenomicAlignments::readGAlignments(file = bam_file,
                                                  param = param)
        
        #' (1) retain reads overlaps with the peaks (must be n:1 mapping, n>=1)

        ol <- .bestOverlaps(query=gap, subject=peaks, ignore.strand=TRUE)
        sel_gap <- gap[queryHits(ol)]
        mcols(sel_gap)$Overlaps.peak.names <- names(peaks)[subjectHits(ol)]
        sel_gap <- sel_gap[njunc(sel_gap) > 0]
        #' (2) retain reads that has junctions
        junc <- junctions(sel_gap)
        tmp <- rep(mcols(sel_gap)$Overlaps.peak.names, times=elementNROWS(junc))
        junc <- unlist(junc)
        mcols(junc)$Overlaps.peak.names <- tmp
        #' (3) retain reads splicing into annotated genes
        ol <- .bestOverlaps_grl(query=junc, subject=tx,
                            ignore.strand=TRUE, minoverlap=20L)
        target_junc <- junc[queryHits(ol)]
        mcols(target_junc)$Overlaps.GeneID <- names(tx)[subjectHits(ol)]
        target_junc
    }, mc.cores=cores, mc.preschedule=FALSE)
    names(list_link_reads) <- NULL
    keep <- elementNROWS(list_link_reads) > 0
    link_reads <- do.call(c, list_link_reads[keep])
    mcols(link_reads)$pair <- paste0(mcols(link_reads)$Overlaps.peak.names, "-",
                                     mcols(link_reads)$Overlaps.GeneID)
    mcols(link_reads)$pair <- factor(mcols(link_reads)$pair)
    link_reads
}



linkPeaksToTx <- function(bam_file, list_which, peaks, tx,
                          mc.cores=1) {
    message(basename(bam_file))
    ## this is for dox samples and find the peak-gene pairs
    list_link_reads <- mclapply(list_which, function(x) {
        flag <- scanBamFlag(isSecondaryAlignment = FALSE)
        param <- ScanBamParam(flag = flag, tag = "XS", which = x)
        gap <- GenomicAlignments::readGAlignments(file = bam_file,
                                                  param = param)
        #' (1) retain reads overlaps with the peaks (must be 1:n mapping, n>=1)
        ##ol <- findOverlaps(gap, peaks, ignore.strand=TRUE)
        ##ol <- removeDupMapping(ol)
        ol <- .bestOverlaps(query=gap, subject=peaks, ignore.strand=TRUE,
                           minoverlap=1L)
        sel_gap <- gap[queryHits(ol)]
        ## need to record with peaks
        mcols(sel_gap)$Overlaps.peak.names <- names(peaks)[subjectHits(ol)]

        #' (2) retain reads link annotated genes
        ##ol <- findOverlaps(sel_gap, tx, ignore.strand=TRUE, minoverlap=20L)
        ##ol <- removeDupMapping(ol)
        ol <- .bestOverlaps_grl(query=sel_gap, subject=tx, ignore.strand=TRUE,
                            minoverlap=20L)
        target_gap <- sel_gap[queryHits(ol)]
        mcols(target_gap)$Overlaps.GeneID <- names(tx)[subjectHits(ol)]
        target_gap
    }, mc.cores=mc.cores, mc.preschedule=FALSE)
    names(list_link_reads) <- NULL
    keep <- elementNROWS(list_link_reads) > 0
    link_reads <- do.call(c, list_link_reads[keep])
    mcols(link_reads)$pair <- paste0(mcols(link_reads)$Overlaps.peak.names, "-",
                                     mcols(link_reads)$Overlaps.GeneID)
    mcols(link_reads)$pair <- factor(mcols(link_reads)$pair)
    link_reads
}


getCommonLinksHits <- function(dox_links, nodox_links) {
    pair_hits <- sapply(dox_links, function(x) {
        table(mcols(x)$pair)
    })
    pair_list <- unique(unlist(sapply(pair_hits, names)))
    hits <- matrix(0, nrow=length(pair_list), ncol=length(dox_links),
                   dimnames=list(pair_list, names(dox_links)))
    for (i in names(pair_hits)) {
        hits[names(pair_hits[[i]]), i] <- pair_hits[[i]]
    }
    dox_link_hits <- hits

    #' get nodox_link_hits
    nodox_link_hits <- .getCommonLinksHits.nodox(dox_link_hits, nodox_links)
    
    list(dox=dox_link_hits, nodox=nodox_link_hits)
}


.getCommonLinksHits.nodox <- function(dox_link_hits, nodox_links) {
    pair_hits <- lapply(nodox_links, function(x) {
        keep <- mcols(x)$pair %in% rownames(dox_link_hits)
        x <- x[keep]
        table(as.character(mcols(x)$pair))
    })
    
    hits <- matrix(0, nrow=nrow(dox_link_hits), ncol=length(nodox_links),
                   dimnames=list(rownames(dox_link_hits), names(nodox_links)))
    for (i in names(pair_hits)) {
        hits[names(pair_hits[[i]]), i] <- pair_hits[[i]]
    }
    hits        
}

linkAnnotation <- function(org, link_hits, splicing_hits, pkgDir, peaks) {
    #' evaluate activation score
    score <- (rowSums(link_hits$dox)+0.5) / (rowSums(link_hits$nodox)+0.5)

    splicing_score <- (rowSums(splicing_hits$dox)+0.5) /
        (rowSums(splicing_hits$nodox)+0.5)

    #' report/annotation - links
    peak_name <- sapply(strsplit(names(score), "-"), "[[", 1)
    gene_id  <- sapply(strsplit(names(score), "-"), "[[", 2)
    gene_name <- select(x=org, keys=gene_id,
                        keytype="ENTREZID", columns=c("SYMBOL"))
    peak_gr <- as(peaks[peak_name], "GRanges")
    df <- DataFrame(Peak.Name=peak_name,
                Peak.Range=as.character(peak_gr),
                Gene.EntrezID=gene_id,
                Gene.Symbol=gene_name[, "SYMBOL"],
                Total.Dox.link.hits=rowSums(link_hits$dox),
                Total.NoDox.link.hits=rowSums(link_hits$nodox),
                Activation.Score=as.numeric(score),
                Peak.ov.LTR=peak_gr$overlap.LTR,
                Peak.repName=peak_gr$repName,
                Peak.repClass=peak_gr$repClass,
                Peak.repFamily=peak_gr$repFamily,
                Splicing.link=names(score) %in% names(splicing_score),
                as(link_hits$dox, "DataFrame"),
                as(link_hits$nodox, "DataFrame"))

    rownames(df) <- names(score)
    df <- df[order(df$Activation.Score, decreasing=TRUE), ]
    df <- df[df$Activation.Score >= 3, ]
    save(df, file=file.path(pkgDir, "data", "PeakActivatedGene.rda"))
    write.csv(df, file=file.path(pkgDir, "inst", "stats", "PeakActivatedGene.csv"))
    df
}

splicingLinkAnnotation <- function(org, link_hits, splicing_hits, pkgDir, peaks) {
    #' evaluate activation score
    score <- (rowSums(link_hits$dox)+0.5) / (rowSums(link_hits$nodox)+0.5)

    splicing_score <- (rowSums(splicing_hits$dox)+0.5) /
        (rowSums(splicing_hits$nodox)+0.5)

    #' report/annotation - links
    peak_name <- sapply(strsplit(names(splicing_score), "-"), "[[", 1)
    gene_id   <- sapply(strsplit(names(splicing_score), "-"), "[[", 2)
    gene_name <- select(x=org, keys=gene_id,
                        keytype="ENTREZID", columns=c("SYMBOL"))
    peak_gr <- as(peaks[peak_name], "GRanges")
    df <- DataFrame(Peak.Name=peak_name,
                Peak.Range=as.character(peak_gr),
                Gene.EntrezID=gene_id,
                Gene.Symbol=gene_name[, "SYMBOL"],
                Total.Dox.link.hits=rowSums(splicing_hits$dox),
                Total.NoDox.link.hits=rowSums(splicing_hits$nodox),
                Activation.Score=as.numeric(splicing_score),
                Peak.ov.LTR=peak_gr$overlap.LTR,
                Peak.repName=peak_gr$repName,
                Peak.repClass=peak_gr$repClass,
                Peak.repFamily=peak_gr$repFamily,
                hits.link=names(splicing_score) %in% names(score),
                as(splicing_hits$dox, "DataFrame"),
                as(splicing_hits$nodox, "DataFrame"))

    rownames(df) <- names(splicing_score)
    df <- df[order(df$Activation.Score, decreasing=TRUE), ]
    df <- df[df$Activation.Score >= 3, ]
    save(df, file=file.path(pkgDir, "data", "PeakActivatedGene_splicing.rda"))
    write.csv(df, file=file.path(pkgDir, "inst", "stats",
                      "PeakActivatedGene_splicing.csv"))
    df
}
