pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.MB135.HDUX4"

library(DESeq2)
library(ggplot2)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#'
#' load data
#' 
library(hg38.MB135.HDUX4)
data(HinH.knownGene)
se <- HinH.knownGene[rowSums(assay(HinH.knownGene)) > 6, ]

#'
#' Sample Distance PCA
#'
dds <- DESeqDataSet(se, design = ~ Treatment)
rlg <- rlog(dds)
data <- plotPCA(rlg, intgroup=c("Treatment"), returnData=TRUE)
data$Sample.Names <- rlg$Sample.Names
percentVar <- round(100 * attr(data, "percentVar"))
png(file.path(pkgDir, "inst", "figures", "SamplePCAplot.png"))
qplot(PC1, PC2, color=Sample.Names, shape=Treatment, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + labs(title="HinH")
dev.off()

pdf(file.path(pkgDir, "inst", "figures", "SamplePCAplot.pdf"))
qplot(PC1, PC2, color=Sample.Names, shape=Treatment, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + labs(title="HinH")
dev.off()  

#'
#' Spearsman correlation test for replica
#'
corMatrix <- corP <- matrix(NA, ncol=ncol(se), nrow=ncol(se),
                            dimnames=list(colnames(se), colnames(se)))
for (i in 1:ncol(se)) {
    for (j in 1:ncol(se)) {
        tmp <- cor.test(x=assay(se)[, i], y=assay(se)[, j],
                        method="spearman", alternative="two.sided")
        corMatrix[i, j] <- tmp$estimate
        corP[i, j] <- tmp$p.value
    }
}
library(pheatmap)
png(file.path(pkgDir, "inst", "figures", "ImageSpearmanCor.png"))
srho <- min(min(corMatrix[c(1, 3, 5), c(1,3,5)]), min(corMatrix[c(2,4,6), c(2,4,6)]))
pheatmap(corMatrix, main=paste0("min(rho): ", 
                           format(srho, digits=3)),
            margins=c(10,10))
dev.off()

#'
#' DEseq: differential analysis and output results
#' Entrez Gene ID;Entrez Gene Symbol;HGNC Gene Symbol (for human);
#' MGI Gene Symbol (for mouse); Gene name (The written out form of the gene symbol.);
#' HomoloGene ID
#' 
HinH.dds <- estimateSizeFactors(dds)
HinH.dds <- DESeq(dds)
save(HinH.dds, file=file.path(pkgDir, "data", "HinH.dds.rda"))

#'
#' report
#' 
library(org.Hs.eg.db)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/miscTools.R")
HinH.report <- genReport.DDS.hg38(HinH.dds)
 

save(HinH.report, file=file.path(pkgDir, "data", "HinH.report.rda"))

convDFelement <- function(DF) {
    ## DF is the DataFrame reort
    chr <- sapply(DF$chr, function(x) as.character(x)[1])
    start <- sapply(DF$start, function(x) paste(as.character(x), collapse=";"))
    end <- sapply(DF$end, function(x) paste(as.character(x), collapse=";"))
    strand <- sapply(DF$strand, function(x) paste(as.character(x), collapse=";"))
    DF$chr <- chr
    DF$start <- start
    DF$end <- end
    DF$strand <- strand
    as.data.frame(DF)
}

df <- convDFelement(HinH.report)
write.csv(df, file=file.path(pkgDir, "inst", "stats", "HinH.report.csv"))

#'
#' Repeats analysis: 229 repeats and 33 corresponding families are enriched
#'
library(DESeq2)
load(file.path(pkgDir, "data", "SE.repName.rda"))
se <- SE.repName[rowSums(assay(SE.repName)) > 6, ]
dds.repName <- DESeqDataSet(se, design = ~ Treatment)
dds.repName <- estimateSizeFactors(dds.repName)
dds.repName <- DESeq(dds.repName)
res <- results(dds.repName)
sig <- which(res$padj < 0.001)
head(res[sig, ])

repFam <- sapply(rowRanges(se), function(x) {
    fam <- factor(x$repFamily)
    paste(levels(fam), collapse=",")
})

unique(repFam[sig])


## have a list first:
tmp <- rep("No", length(repFam))
tmp[sig] <- "Yes"

dat <- data.frame(family=repFam,
                  Sig=tmp)
## take out some rows in which teh family has zero sig
tb.sig <- table(repFam[sig])
tb.all <- table(repFam)
p <- names(tb.all) %in% names(tb.sig)
dontkeep <- names(tb.all[!p])
i <- dat$family %in% dontkeep
dat <- dat[!i, ]

ratio.str <- paste0(as.character(tb.sig), "/", as.character(tb.all[names(tb.sig)]))

png("test.png")
ggplot(dat, aes(family, fill=sig)) + geom_bar() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

png("test2.png")
barplot(tb.sig/tb.all[names(tb.sig)], las=2, cex.axis=0.5)
dev.off()

#'
#' Enrichment testing? universe=1226, selected=299
#' hypergeomatric testing for the repFamily enrichement analysis?
#' 
universe=names(repFam)
selected=names(repFam[sig])
tb.all <- table(repFam)
tb.sig <- table(repFam[sig])
cnt <- rep(0, length(tb.all))
names(cnt) <- names(tb.all)
cnt[names(tb.sig)] <- tb.sig

m <- 299
n <- 997
x <- cnt
k <- as.numeric(tb.all)
prob <- dhyper(x=cnt, m=m, n=n, k=k)
prob[prob < 0.1] ## and cnt > mu to be over-represented

mu <- k*(m/(m+n))
