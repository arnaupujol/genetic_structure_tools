### Create data.frame with reads
library(dplyr)
library(tidyr)
library(haven)
fold <- "results/MDA_R3/"

load(paste0(fold, "demultiplexMarkerSummary.Rdata"))
load(paste0(fold, "HaplotypesTable_rep1_s3.Rdata"))
pheno <- read_dta("data/DB_Carlos_19Oct2019.dta")

## Get common samples
pheno.f <- subset(pheno, !nida %in% c(1661030.9, 1648537.2))
CspHapTable <- filter(tab_rep1_s3$Csp, !is.na(Haplotype))
AmaHapTable <- filter(tab_rep1_s3$Ama, !is.na(Haplotype))
comSamples <- Reduce(intersect, list(pheno.f$sample, 
                                     CspHapTable$SampleName, 
                                     AmaHapTable$SampleName))

CspHapTable.f <- subset(CspHapTable, SampleName %in% comSamples)
AmaHapTable.f <- subset(AmaHapTable, SampleName %in% comSamples)

CspMerge <- CspHapTable.f %>% 
  group_by(SampleName) %>% 
  summarize(CspMergeReads = sum(Reads))
AmaMerge <- AmaHapTable.f %>% 
  group_by(SampleName) %>% 
  summarize(AmaMergeReads = sum(Reads))

readsSummary <- deplexMarker %>%
  select(SampleName, MarkerID, numReadIn, numReadOut) %>% 
  spread(MarkerID, numReadOut) %>% 
  inner_join(AmaMerge, by = "SampleName") %>%
  inner_join(CspMerge, by = "SampleName") %>%
  mutate(InitialReads = numReadIn,
         DemultiplexCspReads = Csp, 
         DemultiplexAmaReads = Ama) %>%
  select(SampleName, ends_with("Reads"))
write.table(readsSummary, file = paste0(fold, "summaryReads.txt"), quote = FALSE, row.names = FALSE)
