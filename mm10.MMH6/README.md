# MMH6

## Description
Triplicate of MMH6 RNA-Seq samples in DOX and NO-DOX.

## Platform
Illumina HS2500? 

## Overall Design (GEO)
We generated a new RNA-seq dataset in mouse myoblasts and compared
them to a published dataset for hDUX4 in human myoblasts. Our MMH
(mDUX homeodomains with hDUX4 c-term) cells were polyclonal. The
dataset (and the dataset we compare to) used a doxycycline-inducible
system and a codon-altered transgene. It reflect biological
triplicates of clonal cell lines induced with doxycycline for
thirty-six hours; doxycycline induction caused expression of mDUX in
the mouse myoblast dataset. The control datasets, prepared in
biological triplicates,  were uninduced cells from the MMH cell
lines. 

## Data Origin
- MMH6: `/shared/ngs/illumina/jwhiddon/161121_D00300_0344_AH5MTNBCXY`
- MinM and Luciferase: `/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX`

## Aligned Bams
- MMH6: `/shared/ngs/illumina/jwhiddon/161121_D00300_0344_AH5MTNBCXY/tophat/bam`
- MinM and Luciferase: `/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX/tophat/bam`

## Genome Build
- mm10

## Software and Versions
- R 3.3/Bioconductor 3.3
- tophat-2.1.0/bowtie2-2.2.3

## Annotation Track
- UCSC KnownGene (Bioconductor:TxDb.Mmusculus.UCSC.mm10.knownGene)

## Analysis Package
`/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.MMH6`


### Directory
1. data:
  - `MMH6.dds`: MMH6 DOX vs MMH6 NODOX
  - `MMH6_luc.dds`: DOX vs Luc NODOX
  - `MinM.dds`: MinM DOX vs MinM NODOX
  -  `MinM_luc.dds`: MinM DOX vs Luc NODOX

2. inst/stats:
  - `MMH6_MinM_lfc.csv`: containing the annotation and log fold change of MMH6-DOX (vs MMH6-NODOX) and MinM-DOX (vs MinM-NODOX).
  - `MMH6_MinM__Luc_lfc.csv`: containing the annotation and log fold change of MMH6 (vs Luc-NODOX) and MinM (vs Luc-NODOX). The lfc between MMH6 and MinM gives better correlation than the previous comparison. 

NOTE that the padj value is based upon the hypothesis that |lfc| > 2. The threshold for padj to call differential expression is 0.1.

## UCSCTrackHub
- UCSC track name: `C2C12_MMH6`
- bigWig files: `/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/mm10_C2C12_MMH6_RNASeq`

## Workflows
1. Preprocess: Alignement using tophat and bowtie2
2. Count Hits: "Intersection-strict" mode
3. DESeq2
  - MMH6 DOX vs MMH6 NODOX
  - MMH6 DOX vs Luc NODOX
  - MinM DOX vs MinM NODOX
  - MinM DOX vs Luc NODOX

 
