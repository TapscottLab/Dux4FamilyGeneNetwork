#!/app/easybuild/software/R/3.3.1-foss-2016b-fh1/bin/Rscript
#./analysis.R
#SBATCH -N1 -n6 -t 2-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

##library(processSeqTools)
## module load R/3.3.1-foss-2016b-fh1; to test the mm10.knownGene track is the one that built in Oct 05, 2015

workers <- 6L
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
source("/fh/fast/tapscott_s/CompBio/R_package/processSeqTools/R/makeSE.R")
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/simpleDESeqReport.R")

#'
#' MMH6 and DESeq DOX vs NODOX
#' 
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
ngsDir <- "/shared/ngs/illumina/jwhiddon/161121_D00300_0344_AH5MTNBCXY"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$",  include.dirs=TRUE,
                       full.names=TRUE)
Treatment <- factor(c(rep("NODOX",3), rep("DOX", 3)), levels=c("NODOX", "DOX"))
res <- makeSEwrapper(bamFiles=bamFiles, pkgDir=pkgDir,
                     OrgDb=org.Mm.eg.db,
                     TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     features=NULL, workers=6L,
                     Treatment=Treatment, title="MMH6",
                     plot.PCA=TRUE,
                     do.DESeq=TRUE, print.report=TRUE)

#' some process to change the names of se and dds
setwd(file.path(pkgDir, "data"))
load("se.rda")
MMH6.knownGene <- se
save(MMH6.knownGene, file="MMH6.knownGene.rda")
load("dds.rda")
MMH6.dds <- dds
save(MMH6.dds, file="MMH6.dds.rda")

system("rm se.rda")
system("rm dds.rda")

#'
#' make Track
#'
library(UCSCTrackTools)
bwDir <- "/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/mm10_C2C12_MMH6_RNASeq"
makeCoverageTracksWrapper(bamFile=bamFiles, bwDir=bwDir, singleEnded=TRUE,
                  NH.weight=FALSE, cores=workers,
                  trackName="C2C12_MMH6_RNASeq",
                  col = c(102, 194, 165), pattern = "\\.bw$",
                  shortLabel = NULL, longLabel ="C2C12_MMH6_RNASeq Coverage")

#'
#' MinM (DOX vs NODOX)
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
ngsDir <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
bamFiles <- bamFiles[grep("C2C12_mDux", bamFiles)]
Treatment <- factor(c(rep("DOX",3), rep("NODOX", 3)), levels=c("NODOX", "DOX"))
res <- makeSEwrapper(bamFiles=bamFiles, pkgDir=pkgDir,
                     OrgDb=org.Mm.eg.db,
                     TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     features=NULL, workers=6L,
                     Treatment=Treatment, title="MinM",
                     plot.PCA=FALSE,
                     do.DESeq=TRUE, print.report=TRUE)
setwd(file.path(pkgDir, "data"))
load("se.rda")
MinM.knownGene <- se
save(MinM.knownGene, file="MinM.knownGene.rda")
load("dds.rda")
MinM.dds <- dds
save(MinM.dds, file="MinM.dds.rda")
system("rm se.rda")
system("rm dds.rda")

#'
#' Luciferase and DESeq (DOX vs NODOX)
#'
source("~/tapscott/R_package/processSeqTools/R/makeSE.R")
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
ngsDir <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
bamFiles <- bamFiles[grep("Luc", bamFiles)]
Treatment <- factor(c(rep("DOX",3), rep("NODOX", 3)), levels=c("NODOX", "DOX"))
res <- makeSEwrapper(bamFiles=bamFiles, pkgDir=pkgDir,
                     OrgDb=org.Mm.eg.db,
                     TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     features=NULL, workers=6L,
                     Treatment=Treatment, title="Luc",
                     plot.PCA=FALSE,
                     do.DESeq=TRUE, print.report=TRUE)
setwd(file.path(pkgDir, "data"))
load("se.rda")
Luciferase.knownGene <- se
save(Luciferase.knownGene, file="Luciferase.knownGene.rda")
load("dds.rda")
Luciferase.dds <- dds
system("rm se.rda")
system("rm dds.rda")
      
#'
#' Compare MinM (dox/no-dox) and MMH6 (dox/no-dox)
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
dataDir <- file.path(pkgDir, "data")
load(file.path(dataDir, "MMH6.dds.rda"))
load(file.path(dataDir, "MinM.dds.rda"))

library(DESeq2)
library(ggplot2)
MMH6.res <- results(MMH6.dds, lfcThreshold=2, alpha=0.1)
MinM.res <- results(MinM.dds, lfcThreshold=2, alpha=0.1)
genes <- intersect(rownames(MMH6.res), rownames(MinM.res))
MMH6.res <- MMH6.res[genes, ]
MinM.res <- MinM.res[genes, ]
df <- data.frame(MMH6.lfc = MMH6.res$log2FoldChange,
                 MinM.lfc = MinM.res$log2FoldChange,
                 MMH6.sig = MMH6.res$padj < 0.1,
                 MinM.sig = MinM.res$padj < 0.1)
df$DE <- "None"
df$DE[df$MMH6.sig  & df$MinM.sig] <- "Both"
df$DE[df$MMH6.sig  & !df$MinM.sig] <- "MMH6 only"
df$DE[!df$MMH6.sig & df$MinM.sig] <- "MinM only"
rownames(df) <- genes
df$sym <- mapIds(org.Mm.eg.db, keys=genes, colum="SYMBOL", keytype="ENTREZID", multiVals="first")
write.csv(df, file=file.path(pkgDir, "inst", "stats", "MMH6_MinM_lfc.csv"))

cor <- cor(df$MMH6.lfc, df$MinM.lfc, method="pearson") ## 0.755
lm <- summary(lm(df$MMH6.lfc ~ df$MinM.lfc)) ## y=0.004532 + 0.864044x; R2=0.6
png(file.path(pkgDir, "inst", "figures", "scatter_MMH6_MinM_lfc.png"))
ggplot(df, aes(x=MMH6.lfc, y=MinM.lfc, color=DE)) + geom_point(size=0.5, alpha=0.5) +
    geom_rug() + annotate("text", x=-3, y=10, label="Pearson Cor: 0.755\nR2: 0.57",
                           hjust=0, vjust=0, color="black", size=4)
dev.off()

png(file.path(pkgDir, "inst", "figures", "scatter_MMH6_MinM_lfc_BW.png"))
ggplot(df, aes(x=MMH6.lfc, y=MinM.lfc)) + geom_point(size=0.5, alpha=0.5) +
    annotate("text", x=-3, y=10, label="Pearson Cor: 0.755\nR2: 0.57",
                           hjust=0, vjust=0, color="black", size=4)
dev.off()
table(df$DE)

#'
#' Compare MMH6 and MinM (Luc as no-dox) (dox/luc no-dox and MinM dox/no-dox) using
#' Luciferase samples as control
#' 
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
setwd(file.path(pkgDir, "data"))
load("MMH6.knownGene.rda")
load("MinM.knownGene.rda")
load("Luciferase.knownGene.rda")
library(DESeq2)
library(ggplot2)

## sample distance
figDir <- file.path(pkgDir, "inst", "figures")
se <- cbind(MMH6.knownGene[, MMH6.knownGene$Treatment=="DOX"],
            MinM.knownGene[, MinM.knownGene$Treatment=="DOX"],
            Luciferase.knownGene[, Luciferase.knownGene$Treatment=="NODOX"])
se$Nickname <- factor(c(rep("MMH6-DOX", 3), rep("MinM-DOX", 3), rep("Luc-NODOX", 3)))
dds <-  DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ Nickname)
rlg <- rlog(dds)
data <- plotPCA(rlg, intgroup=c("Nickname"), returnData=TRUE) 
percentVar <- round(100 * attr(data, "percentVar"))
title <- "MMH-MinM-Luc"
png(file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
qplot(PC1, PC2, color=Nickname, data=data) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
           ylab(paste0("PC2: ",percentVar[2],"% variance")) + labs(title=title)
dev.off()

## DESeq
se <- cbind(MMH6.knownGene[, MMH6.knownGene$Treatment=="DOX"],
            Luciferase.knownGene[, Luciferase.knownGene$Treatment=="NODOX"])
MMH6.dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ Treatment)
MMH6.dds <- DESeq(MMH6.dds)
save(MMH6.dds, file=file.path(pkgDir, "data", "MMH6_Luc.dds"))
     
se <- cbind(MinM.knownGene[, MinM.knownGene$Treatment=="DOX"],
            Luciferase.knownGene[, Luciferase.knownGene$Treatment=="NODOX"])
MinM.dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ Treatment)
MinM.dds <- DESeq(MinM.dds)
save(MinM.dds, file=file.path(pkgDir, "data", "MinM_Luc.dds"))

## comparison
library(DESeq2)
library(ggplot2)
library(org.Mm.eg.db)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6"
figDir <- file.path(pkgDir, "inst", "figures")
load(file.path(pkgDir, "data", "MMH6_Luc.dds"))
load(file.path(pkgDir, "data", "MinM_Luc.dds"))

MMH6.res <- results(MMH6.dds, lfcThreshold=2, alpha=0.1)
MinM.res <- results(MinM.dds, lfcThreshold=2, alpha=0.1)
genes <- intersect(rownames(MMH6.res), rownames(MinM.res))
MMH6.res <- MMH6.res[genes, ]
MinM.res <- MinM.res[genes, ]

df <- data.frame(MMH6.lfc = MMH6.res$log2FoldChange,
                 MinM.lfc = MinM.res$log2FoldChange,
                 MMH6.sig = MMH6.res$padj < 0.1,
                 MinM.sig = MinM.res$padj < 0.1)
df$DE <- "None"
df$DE[df$MMH6.sig  & df$MinM.sig] <- "Both"
df$DE[df$MMH6.sig  & !df$MinM.sig] <- "MMH6 only"
df$DE[!df$MMH6.sig & df$MinM.sig] <- "MinM only"
rownames(df) <- genes
df$sym <- mapIds(org.Mm.eg.db, keys=genes, colum="SYMBOL", keytype="ENTREZID", multiVals="first")
write.csv(df, file=file.path(pkgDir, "inst", "stats", "MMH6_MinM_Luc_lfc.csv"))

cor <- cor(df$MMH6.lfc, df$MinM.lfc, method="pearson") ## 0.78
lm <- summary(lm(df$MMH6.lfc ~ df$MinM.lfc)) ## y=0.129597x + 0.896912; R2=0.6158
png(file.path(pkgDir, "inst", "figures", "scatter_MMH6_MinM_Luc_lfc.png"))
ggplot(df, aes(x=MMH6.lfc, y=MinM.lfc, color=DE)) + geom_point(size=0.5, alpha=0.5) +
    geom_rug() + annotate("text", x=-5, y=10, label="Pearson Cor: 0.784 \nR2=0.6158",
                           hjust=0, vjust=0, color="black", size=4)
dev.off()

png(file.path(pkgDir, "inst", "figures", "scatter_MMH6_MinM__Luc_lfc_BW.png"))
ggplot(df, aes(x=MMH6.lfc, y=MinM.lfc)) + geom_point(size=0.5, alpha=0.5) +
    annotate("text", x=-5, y=10, label="Pearson Cor: 0.784 \nR2=0.6158",
                           hjust=0, vjust=0, color="black", size=4)
dev.off()
