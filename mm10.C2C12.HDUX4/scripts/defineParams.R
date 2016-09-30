pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.HDUX4"
ngsDir <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
bamFiles <- bamFiles[grep("C2C12_HDUX4", bamFiles)]

