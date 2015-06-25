## Evaluation of fastq and bam file qc tools available in bioconductor


############### gc fastq
# seqTools - gcContentMatrix
library(seqTools)
seqFiles <- list.files("~/Desktop/pepr-data/MG002/fastq", 
           pattern = "fastq", full.names = TRUE)
testFq <- fastqq(c(list.files("~/Desktop/pepr-data/MG002/fastq", 
                              pattern = "fastq", full.names = TRUE)))
## issue with memory overflow for PacBio data
pacbioFq <- fastqq(c(seqFiles[41:43]))
plotGCcontent(testFq)
# %GC consistent between datasets

for(i in 1:20){
    plotNucFreq(testFq,i, maxx = 400)
}
# Nucleotide Frequency consistent across read lengths

## test SeqLength
## can use to get seq length values for table
seqLen <- seqLenCount(testFq) 
# use to get min and max - then will need to calculate median
reads <- nReads(testFq)

############### seq length
## maybe qrqc????
##


## Looking into visualization
library(Gviz)
afrom <- 2960000
ato <- 3160000
bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE),stacking = "dense")
alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "gapped.bam"), isPaired = FALSE)
plotTracks(c(bmt, alTrack), from = afrom + 12700,to = afrom + 15200, chromosome = "chr12")
plotTracks(alTrack, from = afrom + 12700,to = afrom + 15200, chromosome = "chr12")
