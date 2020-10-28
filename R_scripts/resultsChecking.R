#'#################################################################################
#'#################################################################################
#' Check obtained results
#'#################################################################################
#'#################################################################################

## Check samples with MOI of 0 ####
outputDir <- "results/MDA_R3/"
load(file = file.path(outputDir, "demultiplexMarkerSummary.Rdata"))

MOI <- read.table(file.path(outputDir, "sample_MOI.txt"), header = TRUE)

summary(subset(deplexMarker, SampleName %in% MOI[is.na(MOI$Ama), "SampleName"] & 
                 MarkerID == "Ama")$numReadOut)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 13.00   20.00   25.00   27.70   35.75   52.00
summary(subset(deplexMarker, SampleName %in% MOI[is.na(MOI$Csp), "SampleName"] & 
                 MarkerID == "Csp")$numReadOut)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 6       7      12      12      15      22
