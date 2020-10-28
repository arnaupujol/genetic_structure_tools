#'#################################################################################
#'#################################################################################
#' Run HaplotypR on dilution samples alone to estimate parameters
#'#################################################################################
#'#################################################################################

### Merge fastq files corresponding to the same sample from different batches (bash)
mkdir results/Dillution/adapted_fastq
out=results/Dillution/adapted_fastq

genes=( csp ama )
for gen in "${genes[@]}"
do
  for f in `ls data/*$gen*` 
  do
    #echo $f
    filename=$(basename $f)
    filename="${filename/R1_001/F}" ## Sustituir valores en variables - como gsub
    filename="${filename/R2_001/R}"
    cp $f $out/$filename
    
  done
done

## Load libraries ####
library(HaplotypR)
library(Biostrings)
library(ShortRead)
library(parallel)


# Set output directory
outputDir <- "results/Dillution/"

## Load markers table ####
markerTab <- read.delim("results/HaplotypR_files/MOI_markerFile.txt", stringsAsFactors = F)

## Create files data.frame ####
dilutionF <- dir(paste0(outputDir, "adapted_fastq"), pattern = "Set.*F", full.names = TRUE)
dilutionR <- dir(paste0(outputDir, "adapted_fastq"), pattern = "Set.*R", full.names = TRUE)
mean(substring(dilutionF, 1, 71) == substring(dilutionR, 1, 71))
# [1] 1 
# File names are correctly mapped

dePlexSample <- data.frame(FileR1 = dilutionF, FileR2 = dilutionR)
parts <- strsplit(dilutionF, split = "_")
genes <- sapply(parts, `[`, 3)
genes <- gsub("fastq/", "", genes)
dePlexSample$SampleID <- sapply(parts, function(x) paste(x[c(4:7)], collapse = "-"))
dePlexSample$SampleID <- paste(genes, dePlexSample$SampleID, sep = "-")
dePlexSample$SampleName <- sapply(parts, `[`, 4)
dePlexSample$SampleName <- paste(genes, dePlexSample$SampleName, sep = "-")

## Add gene to samples data.frame
dePlexMarker <- dePlexSample
dePlexMarker$MarkerID <- ifelse(grepl("CM", genes), "Csp", "Ama")


# Fuse reads ####
### create output subdirectory 
outProcFiles <- file.path(paste0(outputDir, "processedReads"))
dir.create(outProcFiles)

### Set trimming parameters
numNtF <- 170
numNtR <- 120
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

### Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), 
                             substr(refSeq, nchar(refSeq) + 1 - numNtR, nchar(refSeq)), 
                             sep = ""))
names(refSeq) <- markerTab$MarkerID

lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

### Fuse paired read
procReads <- bindAmpliconReads(as.character(dePlexMarker$FileR1), 
                               as.character(dePlexMarker$FileR2), 
                               outProcFiles, 
                               read1Length = numNtF, 
                               read2Length = numNtR)
procReads <- cbind(dePlexMarker[, c("SampleID", "SampleName", "MarkerID")], 
                   procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), 
            sep = "\t", row.names = F)


# Call SNPs ####
### Define function to process SNPs
getSNPsList <- function(marker, minMMrate, minOccGen, postfix, outputDir){
  # Calculate mismatch rate
  files <- as.character(procReads[procReads$MarkerID == marker, "ReadFile"])
  seqErrLst <- calculateMismatchFrequencies(files, 
                                            refSeq[marker], 
                                            method = "pairwiseAlignment",
                                            minCoverage = 100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l) {
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, 
              file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), 
              sep = "\t",
              row.names = F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate = minMMrate, minReplicate = minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start = snp, width = 1))
  }))
  snps <- cbind(Chr = marker, Pos = potSNP, Ref = snpRef, Alt = "N")
  write.table(snps, 
              file = file.path(outputDir, 
                               sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                       minMMrate*100, minOccGen, marker, postfix
                               )
              ), 
              row.names = F, 
              col.names = T, 
              sep = "\t",
              quote = F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, 
                sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                        minMMrate*100, minOccGen, marker, postfix
                )
  ), 
  width = 1500,
  height = 600)
  matplot(seqErr, type = "p", pch = 16, cex = 0.4, col= "#00000088", ylim = c(0, 1),
          ylab = "Mismatch Rate", xlab = "Base Position", main = marker,
          cex.axis = 2, cex.lab = 2)
  abline(v = snps[, "Pos"], lty = 2, col = "grey")
  abline(h = minMMrate, lty = 1, col = "red")
  dev.off()
  return(snps)
}

# Parameters based on SNP calling report
### SNPs Freq: 0.5 
### N samples: 2
snpLstFin <- lapply(markerTab$MarkerID, getSNPsList, 
                   minMMrate = 0.5, minOccGen = 2, postfix = postfix, 
                   outputDir = outputDir)
names(snpLstFin) <- markerTab$MarkerID




# Call haplotype options ####
### Default
minCov <- 3
detectionLimit <- 1e-5
minOccHap <- 6
minCovSample <- 25
options("HaplotypR.devel"=  TRUE) ## Set devel option to get more files

snpLst <- snpLstFin

finalTab <- createFinalHaplotypTable(outputDir = outputDir, 
                                     sampleTable = procReads, 
                                     snpList = snpLst, 
                                     postfix = postfix, 
                                     refSeq = refSeq,
                                     minHaplotypCoverage = minCov, 
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit, 
                                     minSampleCoverage = minCovSample)



## Search haplotypes
CspMap <- read.delim("results/Dillution/HaplotypeOverviewTable_Csp_bind300_275.txt",
                          as.is = TRUE)
CspUID <- read.delim("results/Dillution/contingencyTable_Csp_bind300_275.txt",
                     as.is = TRUE)
cspSeqs <- readDNAStringSet("results/Dillution/frequencyFiles/allSequences_Csp_bind300_275.fasta")
cspHaps <- CspUID %>%
  mutate(total = rowSums(.),
         HaplotypesName = rownames(.)) %>%
  dplyr::select(total, HaplotypesName) %>%
  right_join(CspMap, by = "HaplotypesName") %>%
  filter(FinalHaplotype %in% rownames(CspHap)) %>%
  group_by(FinalHaplotype) %>%
  summarize(maxProp = max(total/sum(total)),
            maxUID = HaplotypesName[which.max(total)])

cspBlast <- cspSeqs[cspHaps$maxUID]
names(cspBlast) <- cspHaps$FinalHaplotype
writeFasta(cspBlast, file = "results/Dillution/Csp_Haplos.fasta")



AmaMap <- read.delim("results/Dillution/HaplotypeOverviewTable_Ama_bind300_275.txt",
                     as.is = TRUE)
AmaUID <- read.delim("results/Dillution/contingencyTable_Ama_bind300_275.txt",
                     as.is = TRUE)
amaSeqs <- readDNAStringSet("results/Dillution/frequencyFiles/allSequences_Ama_bind300_275.fasta")
amaHaps <- AmaUID %>%
  mutate(total = rowSums(.),
         HaplotypesName = rownames(.)) %>%
  dplyr::select(total, HaplotypesName) %>%
  right_join(AmaMap, by = "HaplotypesName") %>%
  filter(FinalHaplotype %in% rownames(AmaHap)) %>%
  group_by(FinalHaplotype) %>%
  summarize(maxProp = max(total/sum(total)),
            maxUID = HaplotypesName[which.max(total)])

amaBlast <- amaSeqs[amaHaps$maxUID]
names(amaBlast) <- amaHaps$FinalHaplotype
writeFasta(amaBlast, file = "results/Dillution/Ama_Haplos.fasta")
