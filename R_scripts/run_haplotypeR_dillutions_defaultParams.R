#'#################################################################################
#'#################################################################################
#' Run HaplotypR on dilution samples alone with predefined parameters
#'#################################################################################
#'#################################################################################

### Merge fastq files corresponding to the same sample from different batches (bash)
mkdir results/DillutionDefault/adapted_fastq
out=results/DillutionDefault/adapted_fastq

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
outputDir <- "results/DillutionDefault/"

## Load markers table ####
markerTab <- read.delim("results/HaplotypR_files/MOI_markerFile_Lerch.txt", stringsAsFactors = F)

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
dePlexSample$BarcodePair <- sapply(parts, function(x) paste(x[c(4:5, 7)], collapse = "-"))


## Demultiplex by marker and remove primers ####
#### create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, trimFilenameExt = "F\\.fastq.gz")

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
save(dePlexMarker, file = file.path(outputDir, "demultiplexMarkerSummary.Rdata"))

## Filter samples with few reads (wrong marker or bad sequencing)
dePlexMarker <- dePlexMarker[dePlexMarker$numReadOut > 25, ]

# Trim reads ####
### create output subdirectory 
trimFiles <- file.path(paste0(outputDir, "trimReads"))
dir.create(trimFiles)

### Set trimming parameters (based on QC report)
numNtF <- 275 ### Reduce in 22 bases - size of primer
numNtR <- 253
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

trimReads <- function (fastqFileR1, fastqFileR2, outputDir, read1Length = NULL,
          read2Length = read1Length, progressReport = message)
{
  if (length(fastqFileR1) != length(fastqFileR2))
    stop("Vector length of fastqFileR1 and fastqFileR2 not identical.")
  tab <- lapply(seq_along(fastqFileR1), function(i) {
    if (!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(fastqFileR1[i]),
                 "and", basename(fastqFileR2[i]), "...")
    progressReport(detail = msg, value = i)
    outputFile1 <- file.path(outputDir, sub("_F\\.fastq.gz",
                                           "_trim_F\\.fastq.gz", basename(fastqFileR1[i])))
    outputFile2 <- file.path(outputDir, sub("_R\\.fastq.gz",
                                            "_trim_R\\.fastq.gz", basename(fastqFileR2[i])))
    f1 <- FastqStreamer(fastqFileR1[i])
    f2 <- FastqStreamer(fastqFileR2[i])
    mode <- "w"
    numReads <- 0
    while (length(sr1 <- yield(f1)) > 0) {
      sr2 <- yield(f2)
      numReads <- numReads + length(sr1)
      if (!is.null(read1Length))
        sr1 <- narrow(sr1, start = 1, width = ifelse(width(sr1) >=
                                                       read1Length, read1Length, NA))
      if (!is.null(read2Length))
        sr2 <- narrow(sr2, start = 1, width = ifelse(width(sr2) >=
                                                       read2Length, read2Length, NA))
      writeFastq(ShortReadQ(sread = sread(sr1),
                            qual = quality(quality(sr1)),
                            id = id(sr1)), file = outputFile1, mode = mode, compress = T)
      writeFastq(ShortReadQ(sread = sread(sr2),
                            qual = quality(quality(sr2)),
                            id = id(sr2)), file = outputFile2, mode = mode, compress = T)
      mode <- "a"
    }
    close(f1)
    close(f2)
    rm(sr1, sr2)
    gc()
    gc()
    if (numReads == 0)
      return(NULL)
    else return(data.frame(numRead = numReads, FileR1 = outputFile1,
                           FileR2 = outputFile2, stringsAsFactors = F))
  })
  tab <- do.call(rbind.data.frame, tab)
  return(tab)
}

trimFiles <- trimReads(as.character(dePlexMarker$FileR1), 
                               as.character(dePlexMarker$FileR2), 
                       trimFiles, 
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
                                 sprintf("SNPs%.0f_occ%i_%s.png", 50, 2, postfix)))



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



finalTab$Ama %>% 
  group_by(SampleName, Haplotype) %>%
  summarize(Reads = sum(Reads)) %>%
  group_by(SampleName) %>%
  mutate(Freq = Reads/sum(Reads),
         Set = substring(SampleName, 8, 10)) %>%
  ggplot(aes(x = SampleName, y = Freq, col = Haplotype)) + 
  geom_point() 


finalTab$Ama %>% 
  group_by(SampleName, Haplotype) %>%
  summarize(Reads = sum(Reads)) %>%
  group_by(SampleName) %>%
  summarize(MOI = sum(grepl("Ama", Haplotype))) %>%
  t()


finalTab$Ama %>% 
  group_by(SampleName, Haplotype) %>%
  summarize(Reads = sum(Reads)) %>%
  group_by(SampleName) %>%
  mutate(Freq = Reads/sum(Reads)) %>%
  filter(grepl("Ama", Haplotype)) %>%
  data.frame()

## Annotate haplotypes (run in ~/data/software/ncbi-blast-2.9.0+-src/c++/ReleaseMT/db)
fold=/scratch/cruiz/32_19_MOI/results/DillutionDefault/
blastn -db nt -query $fold/Csp_HaplotypeSeq_bind170_120.fasta -outfmt   "6 qseqid sseqid stitle pident length mismatch qstart qend evalue bitscore score qcovs qcovus btop" \
-out $fold/Csp_blast.tab -num_threads 40
blastn -db nt -query $fold/Ama_HaplotypeSeq_bind170_120.fasta -outfmt   "6 qseqid sseqid stitle pident length mismatch qstart qend evalue bitscore score qcovs qcovus btop" \
-out $fold/Ama_blast.tab -num_threads 40

## Annotate haplotypes
library(dplyr)
a <- read.delim("res.tab", header = FALSE, as.is = TRUE)
a %>%
  as_tibble() %>%
  group_by(V2) %>%
  summarize(totalScore = sum(V11),
            name = V3[1],
            start1 = min(V7),
            end1 = min(V8),
            length1 = end1 - start1 + 1,
            start2 = max(V7),
            end2 = max(V8),
            length2 = end2 - start2 + 1) %>%
  arrange(desc(totalScore))