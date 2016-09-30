genReport.DDS <- function(dds) {
    res <- results(dds, addMLE=TRUE)
    cnt <- assay(dds)

    symbol  <- mget(rownames(res), org.Mm.egSYMBOL, ifnotfound=NA)
    symbol  <- sapply(symbol, paste, collapse=";")
    name    <- mget(rownames(res), org.Mm.egGeneName, ifnotfound=NA)
    name    <- sapply(name, paste, collapse=";")
    ensembl <- mget(rownames(res), org.Mm.egENSEMBL, ifnotfound=NA)
    ensembl <- sapply(ensembl, paste, collapse=";")

    txs <- range(rowRanges(dds))
    df <- DataFrame(ENTREZ.ID = rownames(dds),
                    symbol=symbol,
                    chr=seqnames(txs),
                    start=start(txs),
                    end=end(txs),
                    strand=strand(txs),
                    res, cnt,
                    Gene.Name=name)
              
}
