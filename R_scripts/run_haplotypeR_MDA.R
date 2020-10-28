#'#################################################################################
#'#################################################################################
#' Run HaplotypR on MDA samples
#'#################################################################################
#'#################################################################################

### Merge fastq files corresponding to the same sample from different batches (bash)
mkdir results/MDA/adapted_fastq
out=results/MDA/adapted_fastq

for f in `ls data/*MDA*` 
do
  #echo $f
  filename=$(basename $f)
  filename="${filename/R1_001/F}" ## Sustituir valores en variables - como gsub
  filename="${filename/R2_001/R}"
  cp $f $out/$filename
done

for f in `ls data/*PM*` 
do
  filename=$(basename $f)
  filename="${filename/R1_001/F}" ## Sustituir valores en variables - como gsub
  filename="${filename/R2_001/R}"
  cp $f $out/$filename
done

## Load libraries ####
### Run everything with default parameters
library(HaplotypR)
library(Biostrings)
library(ShortRead)


# Set output directory
outputDir <- "results/MDA/"

## Load markers table ####
markerTab <- read.delim("results/HaplotypR_files/MOI_markerFile_Lerch.txt", stringsAsFactors = F)

## Create files data.frame ####
mdaF <- dir(paste0(outputDir, "adapted_fastq"), pattern = "F.fastq.gz", full.names = TRUE)
mdaR <- dir(paste0(outputDir, "adapted_fastq"), pattern = "R.fastq.gz", full.names = TRUE)
mean(substring(mdaF, 1, 39) == substring(mdaR, 1, 39))
# [1] 1 
# File names are correctly mapped

dePlexSample <- data.frame(FileR1 = mdaF, FileR2 = mdaR)
parts <- strsplit(mdaF, split = "_")
dePlexSample$SampleID <- sapply(parts, function(x) paste(x[c(4, 8)], collapse = "-"))
dePlexSample$SampleName <- sapply(parts, `[`, 4)
## Add sequencing run to barcode to avoid wrong dePlexMarker tables
dePlexSample$BarcodePair <- sapply(parts, function(x) paste(x[c(5:6, 8)], collapse = "-"))

# Demultiplex by gene ####
# create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, trimFilenameExt = "F\\.fastq.gz")

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)

deplexMarker <-  dePlexMarker[!is.na(dePlexMarker$FileR1),]

# Trim reads ####
### create output subdirectory 
trimFold <- file.path(paste0(outputDir, "trimReads"))
dir.create(trimFold)

### Set trimming parameters (based on QC report)
numNtF <- 275 ### Reduce in 22 bases - size of primer
numNtR <- 253
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

### Load function from run_HaplotypeR_dillutions_defaultParams.R
trimFiles <- trimReads(as.character(dePlexMarker$FileR1), 
                       as.character(dePlexMarker$FileR2), 
                       trimFold, 
                       read1Length = numNtF, 
                       read2Length = numNtR)

# Fuse reads ####
### create output subdirectory 
outProcFiles <- file.path(paste0(outputDir, "processedReads"))
dir.create(outProcFiles)

procReadsMerge <- mergeAmpliconReads(as.character(trimFiles$FileR1), as.character(trimFiles$FileR2), outProcFiles)
procReadsMerge <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReadsMerge)
write.table(procReadsMerge, file.path(outputDir, sprintf("processedReadSummary%s.txt", "_merge")), sep="\t", row.names=F, quote=F)
save(procReadsMerge, file = file.path(outputDir, sprintf("processedReadSummary%s.Rdata", "_merge")))


### Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(refSeq)
names(refSeq) <- markerTab$MarkerID

postfix = "_merge"
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

# Call SNPs ####
### Define function to process SNPs
getSNPsList <- function(marker, minMMrate, minOccGen, procReads, postfix, outputDir){
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

snpLstFin <- lapply(markerTab$MarkerID, getSNPsList, 
                    minMMrate = 0.5, minOccGen = 2, postfix = postfix, 
                    outputDir = outputDir, procReads = procReadsMerge)
names(snpLstFin) <- markerTab$MarkerID
save(snpLstFin, file = file.path(outputDir, 
                                 sprintf("SNPs%.0f_occ%i_%s.Rdata", 50, 2, postfix)))



# Call haplotype options ####
### Default
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2
minCovSample <- 25
snpLst <- snpLstFin

finalTab <- createFinalHaplotypTable(outputDir = outputDir, 
                                     sampleTable = procReadsMerge, 
                                     snpList = snpLst, 
                                     refSeq = refSeq,
                                     postfix = postfix, 
                                     minHaplotypCoverage = minCov, 
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit, 
                                     minSampleCoverage = minCovSample)
save(finalTab, file = file.path(outputDir, "HaplotypesTable.Rdata"))

# Call Haplotype in Run3 samples ####
## Allow haplotype to be only in one sample
## Reduce sensitivity to 1e-5

### Filter procReadsMerge
run3Files <- read.table("data/run3.txt", as.is = TRUE)
parts3 <- strsplit(run3Files$V1, split = "_")
run3Files$BarcodePair <- sapply(parts3, function(x) paste(x[c(4:5, 7)], collapse = "-"))
run3reads <- subset(procReadsMerge, BarcodePair %in% run3Files$BarcodePair)

minCov <- 3
detectionLimit <- 1e-5
minOccHap <- 1
minCovSample <- 25
snpLst <- snpLstFin

tab_rep1_s5 <- createFinalHaplotypTable(outputDir = outputDir, 
                                     sampleTable = run3reads, 
                                     snpList = snpLst, 
                                     refSeq = refSeq,
                                     postfix = postfix, 
                                     minHaplotypCoverage = minCov, 
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit, 
                                     minSampleCoverage = minCovSample)
save(tab_rep1_s5, file = file.path(outputDir, "HaplotypesTable_rep1_s5.Rdata"))


## Preprocess tables ####
### Select haplotypes present in at least 1 sample with a frequency > 1%
### Remove chimeras

### Csp 
CspHapTable <- tab_rep1_s5$Csp

CspFreqTab <- CspHapTable %>% 
  as_tibble() %>%
  filter(!FlagChimera) %>%
  group_by(SampleName, MarkerID) %>%
  mutate(GlobFreq = Reads/sum(Reads), 
         HapFreq = Reads/sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")]),
         Project = ifelse(grepl("PM", SampleName), "PM", "MDA"))

CspHapFreq <- CspFreqTab %>% 
  ungroup() %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  dplyr::select(SampleName, HapFreq, Haplotype) %>%
  spread(SampleName, HapFreq)
CspHapFreq[is.na(CspHapFreq)] <- 0
CspHapsSel <- CspHapFreq$Haplotype[rowMaxs(data.matrix(CspHapFreq[, -1])) > 0.01]

## Ama
AmaHapTable <- tab_rep1_s5$Ama

AmaFreqTab <- AmaHapTable %>% 
  as_tibble() %>%
  filter(!FlagChimera) %>%
  group_by(SampleName, MarkerID) %>%
  mutate(GlobFreq = Reads/sum(Reads), 
         HapFreq = Reads/sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")]),
         Project = ifelse(grepl("PM", SampleName), "PM", "MDA"))

AmaHapFreq <- AmaFreqTab %>% 
  ungroup() %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  dplyr::select(SampleName, HapFreq, Haplotype) %>%
  spread(SampleName, HapFreq)
AmaHapFreq[is.na(AmaHapFreq)] <- 0
AmaHapsSel <- AmaHapFreq$Haplotype[rowMaxs(data.matrix(AmaHapFreq[, -1])) > 0.01]


## Generate MOI table ####
MOI <- rbind(AmaFreqTab, CspFreqTab) %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  filter(Haplotype %in% c(as.character(CspHapsSel), as.character(AmaHapsSel))) %>%
  group_by(SampleName, MarkerID, Project) %>%
  summarize(MOI = length(grep("Ama|Csp", Haplotype))) %>%
  spread(MarkerID, MOI)
write.table(MOI, file = file.path(outputDir, "sample_MOI.txt"), quote = FALSE,
            row.names = FALSE)

## Generate Haplotype Frequency table ####
CspTab <- CspFreqTab %>% 
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  filter(Haplotype %in% CspHapsSel) %>%
  group_by(Haplotype, Project) %>%
  summarize(n = n()) %>%
  spread(Project, n) %>%
  mutate(MDA = ifelse(is.na(MDA), 0, MDA), 
         PM = ifelse(is.na(PM), 0, PM))

write.table(CspTab, file = file.path(outputDir, "Csp_Haps_Freq.txt"), quote = FALSE,
            row.names = FALSE)

AmaTab <- AmaFreqTab %>% 
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  filter(Haplotype %in% AmaHapsSel) %>%
  group_by(Haplotype, Project) %>%
  summarize(n = n()) %>%
  spread(Project, n) %>%
  mutate(MDA = ifelse(is.na(MDA), 0, MDA), 
         PM = ifelse(is.na(PM), 0, PM))


write.table(AmaTab, file = file.path(outputDir, "Ama_Haps_Freq.txt"), quote = FALSE,
            row.names = FALSE)
