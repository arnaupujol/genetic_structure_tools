#'#################################################################################
#'#################################################################################
#' Run HaplotypR on MDA run 3 samples
#'#################################################################################
#'#################################################################################

mkdir results/MDA_R3/
mkdir results/MDA_R3/adapted_fastq
out=results/MDA_R3/adapted_fastq

grep -E 'PM|MDA' data/run3.txt | while read f
do
  filename="${f/R1_001/F}" ## Sustituir valores en variables - como gsub
  filename="${filename/R2_001/R}"
  cp data/$f $out/$filename
done 

## Load libraries ####
### Run everything with default parameters
library(HaplotypR)
library(Biostrings)
library(ShortRead)


# Set output directory
outputDir <- "results/MDA_R3/"

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
dePlexSample$SampleID <- sapply(parts, function(x) paste(x[c(5, 9)], collapse = "-"))
dePlexSample$SampleName <- sapply(parts, `[`, 5)
## Add sequencing run to barcode to avoid wrong dePlexMarker tables
dePlexSample$BarcodePair <- sapply(parts, function(x) paste(x[c(6:7, 9)], collapse = "-"))

### Change wrong sample id
dePlexSample$SampleName[dePlexSample$SampleName == "349"] <- "339"
dePlexSample$SampleID[dePlexSample$SampleID == "349-S295"] <- "339-S295"

# Demultiplex by gene ####
# create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, trimFilenameExt = "F\\.fastq.gz")

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
deplexMarker <-  dePlexMarker[!is.na(dePlexMarker$FileR1),]
save(deplexMarker, file = file.path(outputDir, "demultiplexMarkerSummary.Rdata"))

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



# Custom parameters
## Allow haplotype to be only in one sample
## Reduce sensitivity to 1e-3
minCov <- 3
detectionLimit <- 1e-3
minOccHap <- 1
minCovSample <- 25
snpLst <- snpLstFin

tab_rep1_s3 <- createFinalHaplotypTable(outputDir = outputDir, 
                                     sampleTable = procReadsMerge, 
                                     snpList = snpLst, 
                                     refSeq = refSeq,
                                     postfix = postfix, 
                                     minHaplotypCoverage = minCov, 
                                     minReplicate = minOccHap, 
                                     detectability = detectionLimit, 
                                     minSampleCoverage = minCovSample)
save(tab_rep1_s3, file = file.path(outputDir, "HaplotypesTable_rep1_s3.Rdata"))


## Preprocess tables ####
### Select haplotypes present in at least 1 sample with a frequency > 1%
### Remove chimeras
### Select samples with MOI in both markers

library(dplyr)
library(tidyr)
library(vegan)
library(matrixStats)
library(haven)


load("results/MDA_R3/HaplotypesTable_rep1_s3.Rdata")
pheno <- read_dta("data/DB_Carlos_19Oct2019.dta")

## Process pheno
pheno$study[pheno$study == "MD1"] <- "MDA1"
pheno$study[pheno$study == "MDA1"] <- "MDA"
pheno$study[pheno$study == "XMAG17"] <- "PM"

pheno <- pheno %>%
  mutate(area = as_factor(area),
         area = droplevels(area),
         gender = as_factor(gender),
         fevertemp = as_factor(fevertemp),
         fevertemp = droplevels(fevertemp),
         rdt = as_factor(rdt))

## Remove wrong samples
pheno.f <- subset(pheno, !nida %in% c(1661030.9, 1648537.2))


## Get common samples
CspHapTable <- filter(tab_rep1_s3$Csp, !is.na(Haplotype))
AmaHapTable <- filter(tab_rep1_s3$Ama, !is.na(Haplotype))

comSamples <- Reduce(intersect, list(pheno.f$sample, 
                                     CspHapTable$SampleName, 
                                     AmaHapTable$SampleName))


pheno.f <- subset(pheno.f, sample %in% comSamples)
CspHapTable.f <- subset(CspHapTable, SampleName %in% comSamples)
AmaHapTable.f <- subset(AmaHapTable, SampleName %in% comSamples)


### Csp 
CspFreqTab <- CspHapTable.f %>% 
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

## Select haplotypes present in at least 1 sample with a frequency > 1%
CspHapsSel <- CspHapFreq$Haplotype[rowMaxs(data.matrix(CspHapFreq[, -1])) > 0.01]

write.table(CspHapFreq[CspHapFreq$Haplotype %in% CspHapsSel, ], 
            file = file.path(outputDir, "Csp_Freq_Table_filtered.txt"), quote = FALSE,
            row.names = FALSE)

## Ama
AmaFreqTab <- AmaHapTable.f %>% 
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

## Select haplotypes present in at least 1 sample with a frequency > 1%
AmaHapsSel <- AmaHapFreq$Haplotype[rowMaxs(data.matrix(AmaHapFreq[, -1])) > 0.01]


write.table(AmaHapFreq[AmaHapFreq$Haplotype %in% AmaHapsSel, ], 
            file = file.path(outputDir, "Ama_Freq_Table_filtered.txt"), quote = FALSE,
            row.names = FALSE)



## Add MOI and Shannon to pheno table ####
alphaCsp <- CspHapFreq %>%
  filter(Haplotype %in% CspHapsSel) %>%
  data.matrix() %>%
  t() %>%
  diversity()

alphaAma <- AmaHapFreq %>%
  filter(Haplotype %in% AmaHapsSel) %>%
  data.matrix() %>%
  t() %>%
  diversity()

pheno.f2 <- rbind(AmaFreqTab, CspFreqTab) %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  filter(Haplotype %in% c(as.character(CspHapsSel), as.character(AmaHapsSel))) %>%
  group_by(SampleName, MarkerID) %>%
  summarize(MOI = length(grep("Ama|Csp", Haplotype))) %>%
  spread(MarkerID, MOI) %>%
  mutate(MOI_Ama = Ama, MOI_Csp = Csp) %>%
  dplyr::select(MOI_Ama, MOI_Csp, SampleName) %>%
  right_join(pheno.f, by = c("SampleName" = "sample")) %>%
  mutate(Shannon_Ama = alphaAma[SampleName],
         Shannon_Csp = alphaCsp[SampleName],
         maxFreq_Csp = colMaxs(data.matrix(CspHapFreq[, SampleName])),
         maxFreq_Ama = colMaxs(data.matrix(AmaHapFreq[, SampleName])))

write.table(pheno.f2, file = file.path(outputDir, "sample_Variables.txt"), quote = FALSE,
            row.names = FALSE, sep = "\t")

### Compute genetic relatedness ####
pre <- subset(pheno.f, study == "MDA")$sample
post <- subset(pheno.f, study == "PM")$sample

CspFreqs <- data.matrix(CspHapFreq[CspHapFreq$Haplotype %in% CspHapsSel, -1])
CspFreqs[CspFreqs > 0] <- 1
CspGenRel <- sapply(seq_len(ncol(CspFreqs)), function(i) {
  apply(CspFreqs, 2, function(x) any(CspFreqs[, i] & x))
})
colnames(CspGenRel) <- rownames(CspGenRel)
diag(CspGenRel) <- NA
CspRel <- c(colMeans(CspGenRel[pre, pre], na.rm = TRUE), colMeans(CspGenRel[post, post], na.rm = TRUE))

AmaFreqs <- data.matrix(AmaHapFreq[AmaHapFreq$Haplotype %in% AmaHapsSel, -1])
AmaFreqs[AmaFreqs > 0] <- 1
AmaGenRel <- sapply(seq_len(ncol(AmaFreqs)), function(i) {
  apply(AmaFreqs, 2, function(x) any(AmaFreqs[, i] & x))
})
colnames(AmaGenRel) <- rownames(AmaGenRel)
diag(AmaGenRel) <- NA
AmaRel <- c(colMeans(AmaGenRel[pre, pre], na.rm = TRUE), colMeans(AmaGenRel[post, post], na.rm = TRUE))

pheno.f3 <- pheno.f2 %>%
  mutate(Ama_GenRel = AmaRel[SampleName],
         Csp_GenRel = CspRel[SampleName])

write.table(pheno.f3, file = file.path(outputDir, "sample_Variables.txt"), quote = FALSE,
            row.names = FALSE, sep = "\t")


wilcox.test(colMeans(CspGenRel[pre, pre], na.rm = TRUE), colMeans(CspGenRel[post, post], na.rm = TRUE))

wilcox.test(colMeans(AmaGenRel[pre, pre], na.rm = TRUE), colMeans(AmaGenRel[post, post], na.rm = TRUE))





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


## Generate QC reads table ####
CspHapTable.f %>%
  group_by(SampleName) %>%
  summarize(cov = sum(Reads),
            Hap = sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")])/
              sum(Reads)) %>%
  ungroup() %>%
  dplyr::select(-SampleName) %>%
  summary()
# cov             Hap
# Min.   :  338   Min.   :0.8133
# 1st Qu.: 3872   1st Qu.:0.9704
# Median : 6675   Median :0.9854
# Mean   : 7965   Mean   :0.9737
# 3rd Qu.:10094   3rd Qu.:0.9899
# Max.   :35147   Max.   :0.9958


AmaHapTable.f %>%
  group_by(SampleName) %>%
  summarize(cov = sum(Reads),
            Hap = sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")])/
              sum(Reads)) %>%
  ungroup() %>%
  dplyr::select(-SampleName) %>%
  summary()
# cov             Hap
# Min.   :  307   Min.   :0.7639
# 1st Qu.: 2550   1st Qu.:0.9376
# Median : 4655   Median :0.9767
# Mean   : 6221   Mean   :0.9489
# 3rd Qu.: 8320   3rd Qu.:0.9816
# Max.   :32025   Max.   :0.9881

