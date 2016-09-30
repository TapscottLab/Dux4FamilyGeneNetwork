## This script read the bam files and create a SummarizedExperiments instance containing
## the genomic position and count hits.

library(GenomicFeatures)
library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicAlignments)
library(BiocParallel)
library(ShortRead)
BiocParallel::MulticoreParam()
BiocParallel::registered()

source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4/inst/scripts/defineParams.R")
pkgDir
bamFiles 
###
## Make SummarizedExprients with count reads
###    
genes.mm10 <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene")
HinM.knownGene <- summarizeOverlaps(features=genes.mm10,
                                    reads=BamFileList(bamFiles),
                                    mode="IntersectionStrict",
                                    inter.feature=TRUE,
                                    singleEnd=TRUE,
                                    ignore.strand=TRUE) 

sampleInfo <- data.frame(sample_name = sub(".bam", "", names(bamFiles)),
                         file_bam = path(bamFiles),
                         nick_name = "HinM",
                         Treatment=factor(rep(c("DOX", "NODOX"), each=3),
                             levels=c("NODOX", "DOX")),
                         stringsAsFactors=FALSE)

si <- SGSeq::getBamInfo(sampleInfo, core=6)
save(si, file=file.path(pkgDir, "data", "si.rda"))
colData(HinM.knownGene) <- as(si, "DataFrame")

## FKPM
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/R_scripts/getScaledCounts.R")
assays(HinM.knownGene)$FPKM <- getScaledCountsPerTx(HinM.knownGene)

save(HinM.knownGene, file=file.path(pkgDir, "data", "HinM.knownGene.rda"))

## print out FKPM for Jenn
id <- names(rowRanges(HinM.knownGene))
sym <- select(org.Mm.eg.db, keys=id, keytype="ENTREZID", columns="SYMBOL")
output <- cbind(sym, as.data.frame(assays(HinM.knownGene)$FPKM, "data.frame"))
write.csv(output, file=file.path(pkgDir, "inst", "extdata", "SE.FPKM.csv"))

#'
#' make junction track
#' 

##library(UCSCTrackTools)
source("/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4/inst/scripts/defineParams.R")
source("/fh/fast/tapscott_s/CompBio/R_package/UCSCTrackTools/R/makeJunctionBed.R")
bedDir <- file.path(pkgDir, "inst", "bed")
bbDir <- file.path(pkgDir, "inst", "bigbed")
seqlev <- paste0("chr", c(1:19, "M", "X", "Y"))
makeJunctionBed(bamFiles=bamFiles[1], cores=6, seqlev=seqlev,
                genome="mm10",
                ignore.strand=TRUE, outdir=bedDir, verbose=TRUE)

## Convert bed to bigbed files: create bigbed files with extension .bb
size_file <- "/fh/fast/tapscott_s/CompBio/mm10/mm10.chrom.sizes"
bed_files <- list.files(bedDir, pattern="\\.bed$", full.names=TRUE)
convertBedToBigBed(bed_files, chrom_sizefile=size_file, outdir=bbDir, cores=6)

## move bb files to juncUCSCDir
library(UCSCTrackTools)
bb_files <- list.files(bbDir, pattern="\\.bb$", full.names=TRUE)
juncUCSCDir <- "/home/tapscott/ucsc/junctions/mm10.C2C12.HDUX4"
cmt <- sprintf("cp %s %s", bb_files, juncUCSCDir)
lapply(cmt, system)

## create hub tracklines
createHubTrackLine.junc.bed("mm10.C2C12.HDUX4", hubName="C2C12_HDUX4_Junc")
