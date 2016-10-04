simpleDESeqReport <- function(dds, OrgDb=NULL, thres.padj=NULL, filename) {
    ### dds is the returned values from DESeq
    require(DESeq2)

    if (is.null(OrgDb)) 
        require(org.Hs.eg.db)
    
    if (is.null(thres.padj))
        df <- results(dds)
    else
        df <- subset(results(dds), padj < thres.padj)

    df <- df[order(df$padj, decreasing=FALSE), ]
    ## annotation (EntrzID, ENSEMBL and GENENAME)
    anno <- select(x=OrgDb, keys=rownames(df),
                   keytype="ENTREZID", columns=c("ENSEMBL", "SYMBOL", "GENENAME"))
    anno <- anno[!duplicated(anno$ENTREZID), ]
    rownames(anno) <- anno$ENTREZID
    anno <- anno[rownames(df), ]

    cnt <- counts(dds[rownames(df)])
    output <- cbind(anno, as.data.frame(df), as.data.frame(cnt))
    write.csv(output, file=filename)
    output 
}
