# 1.1.2_Blood_cell_IDOL_estimation


rm(list = ls(all = TRUE))
options(stringsAsFactors=FALSE)
library(minfi)
library(FlowSorted.Blood.EPIC)   # this one should be loaded for EPIC array samples
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kmanifest)  
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
FlowSorted.Blood.EPIC <- libraryDataGet('FlowSorted.Blood.EPIC')
data ("IDOLOptimizedCpGs.compTable")
head(IDOLOptimizedCpGs.compTable)
library(preprocessCore)
# library(genefilter) 	
# library(quadprog)
sessionInfo()

load("/home/avavanasselt/data_working_dir/data/RGset_raw_OLRM.RData")
sampsheet <- read.csv("/home/avavanasselt/data_working_dir/sampsheets/sampsheet_OLRM_2023.csv")
rownames(sampsheet) <- sampsheet$Pool_ID
blood_samps <- sampsheet[sampsheet$Source=="Whole Blood",]
buccal_samps <- sampsheet[sampsheet$Source=="Buccal",]

# RG <- filter(rownames(RG) %in% blood_samps$Pool_ID)
# 
# 
# head(RG)
# str(RG)
#dim(RGset)

# Extract only good quality samples
RGset <- RG[,rownames(blood_samps)]
dim(RGset)

# function estimateCellCounts2.R with the match.arg() line commented out (which produced the error 'arg' should be one of "IlluminaHumanMethylationEPIC")
# source("estimateCellCounts2.R")

annotation(RGset)
annotation(RGset)[which(names(annotation(RGset)) == "array")]
sub("IlluminaHumanMethylation", "", annotation(RGset)[which(names(annotation(RGset)) == "array")])

cellCounts <- estimateCellCounts2(
   rgSet = RGset,
   compositeCellType = "Blood",
   processMethod = "preprocessFunnorm",
   probeSelect = "IDOL",
   cellTypes = c("Bcell", "CD4T", "CD8T", "Neu", "Mono", "NK"),
   referencePlatform = "IlluminaHumanMethylationEPIC",
   referenceset = "FlowSorted.Blood.EPIC",
   IDOLOptimizedCpGs = IDOLOptimizedCpGs.compTable,
   meanPlot = TRUE,
   returnAll = FALSE
 )


save(cellCounts,file="/data/avavanasselt/results/Blood_EPIC_cellcountsidol.RData") 
rm(list = ls(all = TRUE))
gc()


