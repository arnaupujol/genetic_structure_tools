#'#################################################################################
#'#################################################################################
#' Define HaplotypR parameters
#' Process dilution samples with HaplotypeR to define parameters
#'#################################################################################
#'#################################################################################

### Use results from DillutionDefault
mkdir results/HaplotypR_config

## Load libraries ####
library(HaplotypR)
library(Biostrings)
library(ShortRead)

# Set output directory
outputDir <- "results/HaplotypR_config/"

## Load markers table ####
markerTab <- read.delim("results/HaplotypR_files/MOI_markerFile.txt", stringsAsFactors = F)

refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(refSeq)
names(refSeq) <- markerTab$MarkerID

postfix = "_merge"
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})



### Test haplotype calling

load(file = "results/DillutionDefault/processedReadSummary_merge.Rdata")
load("results/DillutionDefault/SNPs50_occ2__merge.png")


### Lower detection limit 
### Use Run3
### Filter procReadsMerge
run3Files <- read.table("data/run3.txt", as.is = TRUE)
parts3 <- strsplit(run3Files$V1, split = "_")
run3Files$BarcodePair <- sapply(parts3, function(x) paste(x[c(3:4, 6)], collapse = "-"))
run3reads <- subset(procReadsMerge, BarcodePair %in% run3Files$BarcodePair)


minCov <- 3
detectionLimit <- 1e-5
minOccHap <- 1
minCovSample <- 25
snpLst <- snpLstFin

options("HaplotypR.devel"=  TRUE) ## Set devel option to get more files

tab_rep1_e5 <- createFinalHaplotypTable(outputDir = outputDir, 
                                     sampleTable = run3reads, 
                                     snpList = snpLstFin, 
                                     postfix = postfix, 
                                     refSeq = refSeq,
                                     minHaplotypCoverage = minCov, 
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit, 
                                     minSampleCoverage = minCovSample)
save(tab_rep1_e5, file = file.path(outputDir, "HaplotypesTable_rep1_e5.Rdata"))
