#'#################################################################################
#'#################################################################################
#' Visualize haplotypes from MOI
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(GenomicAlignments)
library(dplyr)
library(tidyr)

# Convert haplotypes and frequencies to tabular format ####
## Consider haplotype frequency: add one haplotype per 10% frequency in any sample
## Consider presence/absence
## Define functions

makeMarkerTable <- function(freqTab, HapSeq, snps, counts = TRUE){
  
  newFreq <- Reduce(rbind, 
                    lapply(seq_len(nrow(freqTab)), function(x){
                      convertHapsToRows(freqTab[x, ], HapSeq = HapSeq, snps = snps)
                    }))
  if (counts){
    res <- newFreq[rep(seq_len(nrow(newFreq)), newFreq$Freq), ]
  } else {
    res <- newFreq
  }
  res
}

convertHapsToRows <- function(row, HapSeq, snps){
  hap <- HapSeq[as.character(row$Haplotype)]
  genos <- lapply(snps, function(x) paste(rep(substring(hap, x, x), 2), collapse = "/"))
  res <- data.frame(row, genos)
}

load("results/MDA_R3/HaplotypesTable_rep1_s5.Rdata")

### Csp ####
CspHapTable <- tab_rep1_s5$Csp

CspHaps <- readDNAStringSet("results/MDA_R3/Csp_HaplotypeSeq_merge.fasta")
CspSnps <- read.table("results/MDA_R3/potentialSNPlist_rate50_occ2_Csp_merge.txt", header = TRUE)$Pos
names(CspSnps) <- CspSnps

CspFreqs <- CspHapTable %>% 
  as_tibble() %>%
  filter(!FlagChimera) %>%
  filter(!is.na(Haplotype) & !Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")) %>%
  group_by(SampleName) %>%
  mutate(HapFreq = Reads/sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")]),
         Freq = round(Reads/sum(Reads)*10),
         Project = ifelse(grepl("PM", SampleName), "PM", "MDA"))

CspHapFreq <- CspFreqs %>% 
  ungroup() %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  dplyr::select(SampleName, HapFreq, Haplotype) %>%
  spread(SampleName, HapFreq)
CspHapFreq[is.na(CspHapFreq)] <- 0

## Select haplotypes present in at least 1 sample with a frequency > 1%
CspHapsSel <- CspHapFreq$Haplotype[rowMaxs(data.matrix(CspHapFreq[, -1])) > 0.01]
CspFreqs <- filter(CspFreqs, Haplotype %in% CspHapsSel)

CspTabPres <- makeMarkerTable(CspFreqs, CspHaps, CspSnps, counts = FALSE)[ , -c(1:3, 5:8)]

CspTabFreq <- CspFreqs %>%
  group_by(Project, Haplotype) %>%
  mutate(SumFreq = sum(HapFreq)) %>%
  group_by(Project) %>%
  mutate(Freq = round(SumFreq/length(unique(SampleName))*1000)) %>%
  select(Haplotype, Freq, Project) %>%
  distinct()

CspTabCounts <- makeMarkerTable(CspTabFreq, CspHaps, CspSnps)[ , -2]

## Select top haplotype
CspTabTop <- CspFreqs %>% 
  group_by(SampleName) %>% 
  top_n(1, HapFreq) %>%
  select(SampleName, Haplotype, Project) %>%
  makeMarkerTable(., CspHaps, CspSnps, counts = FALSE)

save(CspTabCounts, CspTabPres, CspTabTop, file = "results/MDA_R3/CspTabs.RData")

### Ama ####
AmaHapTable <- tab_rep1_s5$Ama

AmaHaps <- readDNAStringSet("results/MDA_R3/Ama_HaplotypeSeq_merge.fasta")
AmaSnps <- read.table("results/MDA_R3/potentialSNPlist_rate50_occ2_Ama_merge.txt", header = TRUE)$Pos
names(AmaSnps) <- AmaSnps

AmaFreqs <- AmaHapTable %>% 
  as_tibble() %>%
  filter(!FlagChimera) %>%
  filter(!is.na(Haplotype) & !Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")) %>%
  group_by(SampleName) %>%
  mutate(HapFreq = Reads/sum(Reads[!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton")]),
         Freq = round(Reads/sum(Reads)*10),
         Project = ifelse(grepl("PM", SampleName), "PM", "MDA"))

AmaHapFreq <- AmaFreqs %>% 
  ungroup() %>%
  filter(!Haplotype %in% c("Noise", "Indels", "Chimera", "Singelton") &
           !is.na(Haplotype)) %>%
  dplyr::select(SampleName, HapFreq, Haplotype) %>%
  spread(SampleName, HapFreq)
AmaHapFreq[is.na(AmaHapFreq)] <- 0

## Select haplotypes present in at least 1 sample with a frequency > 1%
AmaHapsSel <- AmaHapFreq$Haplotype[rowMaxs(data.matrix(AmaHapFreq[, -1])) > 0.01]
AmaFreqs <- filter(AmaFreqs, Haplotype %in% AmaHapsSel)

AmaTabPres <- makeMarkerTable(AmaFreqs, AmaHaps, AmaSnps, counts = FALSE)[ , -c(1:3, 5:8)]


AmaTabFreq <- AmaFreqs %>%
  group_by(Project, Haplotype) %>%
  mutate(SumFreq = sum(HapFreq)) %>%
  group_by(Project) %>%
  mutate(Freq = round(SumFreq/length(unique(SampleName))*1000)) %>%
  select(Haplotype, Freq, Project) %>%
  distinct()

AmaTabCounts <- makeMarkerTable(AmaTabFreq, AmaHaps, AmaSnps)[ , -2]


## Select top haplotype
AmaTabTop <- AmaFreqs %>% 
  group_by(SampleName) %>% 
  top_n(1, HapFreq) %>%
  select(SampleName, Haplotype, Project) %>%
  makeMarkerTable(., AmaHaps, AmaSnps, counts = FALSE)

save(AmaTabCounts, AmaTabPres, AmaTabTop, file = "results/MDA_R3/AmaTabs.RData")


## Local ####
library(poppr)

load("32_19_MOI/AmaTabs.RData")
load("32_19_MOI/CspTabs.RData")

getMSN <- function(tab){
  gen <- as.genclone(df2genind(tab[, -c(1:2)], sep = "/", ploidy = 1, pop = tab$Project))
  mll.custom(gen) <- as.character(tab$Haplotype)
  gdist <- diss.dist(gen)
  gmsn <- poppr.msn(gen, gdist, showplot = FALSE)
  list(gen = gen, msn = gmsn)
}

# Counts ####
## Csp 
CspCmsn <- getMSN(CspTabCounts)
plot_poppr_msn(CspCmsn$gen, CspCmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1, size.leg = FALSE)

## Ama
AmaCmsn <- getMSN(AmaTabCounts)
plot_poppr_msn(AmaCmsn$gen, AmaCmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1, size.leg = FALSE)


## Presence ####
CspPmsn <- getMSN(CspTabPres)
plot_poppr_msn(CspPmsn$gen, CspPmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1.2, size.leg = FALSE)

AmaPmsn <- getMSN(AmaTabPres)
plot_poppr_msn(AmaPmsn$gen, AmaPmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1.2, size.leg = FALSE)

## Top haplotype ####
CspTmsn <- getMSN(CspTabTop[, -1])
plot_poppr_msn(CspTmsn$gen, CspTmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1.2, size.leg = FALSE)

AmaTmsn <- getMSN(AmaTabTop[, -1])
plot_poppr_msn(AmaTmsn$gen, AmaTmsn$msn, mlg = TRUE, palette = cm.colors, nodescale = 1.2, size.leg = FALSE)


