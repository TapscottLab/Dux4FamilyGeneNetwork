source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4/inst/scripts/defineParams.R")

library(DESeq2)
library(ggplot2)
library(GenomicFeatures)

#'
#' load data
#' 
library(mm10.C2C12.HDUX4)
data(HinM.knownGene)
se <- HinM.knownGene[rowSums(assay(HinM.knownGene)) > 6, ]

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
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

pdf(file.path(pkgDir, "inst", "figures", "SamplePCAplot.pdf"))
qplot(PC1, PC2, color=Sample.Names, shape=Treatment, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()  

#'
#' DEseq: differential analysis and output results
#' Entrez Gene ID;Entrez Gene Symbol;HGNC Gene Symbol (for human);
#' MGI Gene Symbol (for mouse); Gene name (The written out form of the gene symbol.);
#' HomoloGene ID
#' 
dds <- DESeqDataSet(se, design = ~ Treatment)
HinM.dds <- estimateSizeFactors(dds)
HinM.dds <- DESeq(dds)
save(HinM.dds, file=file.path(pkgDir, "data", "HinM.dds.rda"))

## report
library(org.Mm.eg.db)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/miscTools.R")
HinM.report <- genReport.DDS(HinM.dds)

## Adding luciferase log2FC/padj
library(mm10.C2C12.luciferase)
data(luciferase.dds)
luci.res <- results(luciferase.dds)
features <- intersect(rownames(luci.res), rownames(HinM.report))
HinM.report$luciferase.log2FC <- NA
HinM.report$luciferase.padj <- NA
HinM.report[features, "luciferase.log2FC"] <- luci.res[features, "log2FoldChange"]
HinM.report[features, "luciferase.padj"] <- luci.res[features, "padj"]
save(HinM.report, file=file.path(pkgDir, "data", "HinM.report.rda"))


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

tmp <- convDFelement(HinM.report)
write.csv(tmp, file=file.path(pkgDir, "inst", "stats", "HinM.report.csv"))
