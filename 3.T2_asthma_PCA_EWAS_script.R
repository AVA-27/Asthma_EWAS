# Asthma PCA Script

rm(list = ls())

# Load packages
library(rospca)
library(DNAmArray)
library(gee)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(bacon)

# Load methylation data and sample/individual information
load("/data/avavanasselt/data/betas_final.Rdata")
load("/data/avavanasselt/sampsheets/sampsheet_OLRM_2023.RData")
load("/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/pca_scores_final.RData")
smoke_data <- read.csv("/data/avavanasselt/sampsheets/Asthma_epigenetics_Austin_SmokingdataBiobank_30012024.csv", header = TRUE)

# tp1 select 
tp1blood <- sampsheet[rownames(pca_scores),]
table(rownames(tp1blood)==rownames(pca_scores))
rownames(pca_scores) <- tp1blood$FISNumber

# Select TP2 samples
tp2blood <- sampsheet[which(sampsheet$TP2_sample == 1), ]
tp2blood <- tp2blood[tp2blood$Source=="Whole Blood",]
tp2blood <- tp2blood[which(tp2blood$FISNumber %in% rownames(pca_scores)),]
tp2blood <- tp2blood[order(tp2blood$FamilyNumber_los), ]
dim(tp2blood)
# [1] 182 282

#Match sample with correct age at collection
tp2blood <- tp2blood %>%
  mutate(tp2_age = if_else(NTRProject == "BIOBANK1", bioage, b2age))

tp2blood <- merge(tp2blood, smoke_data, by = c("FISNumber", "NTRProject"))



tp2blood <- tp2blood %>%
  mutate(tp2_smoke = if_else(NTRProject == "BIOBANK1", smokstat_bb1, smokstat_bb2))

tp2blood$tp2_smoke

dim(tp2blood)
# [1] 182 298


### Update smoking status ###
# create new column with updated smoking status (0 = never smoke, 1 = previous smoker, 2 = current smoker)
tp2blood <- tp2blood %>%
  mutate(smoke_status = if_else(tp2_smoke == 2, 1, if_else(tp2_smoke == 1, 2, 0)))

table(tp2blood$smoke_status)
#  0  1  2 
#  90 59 33 
# Subset betas based on TP1 samples
tp2blood <- tp2blood[order(tp2blood$FamilyNumber_los), ]
rownames(tp2blood) <- tp2blood$Pool_ID
betas <- betas[, rownames(tp2blood)]

# Transpose betas and align all data
betas <- t(betas)

table(rownames(tp2blood)==rownames(betas))

# subselect PCA data for TP2
pca_tp2 <- pca_scores[which(rownames(pca_scores) %in% tp2blood$FISNumber),]
dim(pca_tp2)
# [1] 182  10
rownames(tp2blood) <- tp2blood$FISNumber
pca_tp2 <- pca_tp2[rownames(tp2blood),]

table(rownames(pca_tp2)==tp2blood$FISNumber)

rownames(tp2blood) <- tp2blood$Pool_ID
table(rownames(tp2blood)==rownames(betas))


# Run EWAS using same PCA data from T1 with new methylation data for T2
npcas <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8_out", "PC9_out", "PC10_out")


load("/home/avavanasselt/data_working_dir/annotation_files/anno_epic_genomicfeatures.RData")

# modify medication covariate due to lack of variability in responses
tp2blood <- tp2blood %>%
  mutate(gen_meds = if_else(med_lungs == 1 | betamimetica1 == 1 | inh_corticosteroid1 == 1, 1, 0))

table(tp2blood$gen_meds)

# Create a data frame for analysis variables
data <- data.frame(
  familynumber = as.numeric(tp2blood$FamilyNumber_los),
  CpGi = rep(NA, nrow(betas)),
  peak_flow = as.numeric(tp2blood$peak_flow),
  sex = as.numeric(tp2blood$sex_los),
  age = as.numeric(tp2blood$tp2_age),
  NK_Perc = as.numeric(tp2blood$NK),
  Mono_Perc = as.numeric(tp2blood$Mono),
  Bcell_Perc = as.numeric(tp2blood$Bcell),
  CD4T_Perc = as.numeric(tp2blood$CD4T),
  CD8T_Perc = as.numeric(tp2blood$CD8T),
  Array_rownum = as.numeric(tp2blood$Array_rownum),
  Sample_Plate = as.factor(tp2blood$Box_plate),
  height = as.numeric(tp2blood$height),
  weight = as.numeric(tp2blood$weight),
  smoke = as.numeric(tp2blood$smoke_status),
  gen_meds = as.numeric(tp2blood$gen_meds)
  
)
rownames(data) <- rownames(tp2blood)

# Make dataframe to record run summary data for each EWAS
run_summary = data.frame()

# Iterate through biomarkers and run EWAS
for (biomarker in npcas) {
  cat("Running EWAS for biomarker:", biomarker, "\n")
  data$tester <- pca_tp2[,biomarker]
  mydata <- data[complete.cases(data[, "tester"]), ]
  samps <- tp2blood[rownames(mydata),]
  samps <- samps[order(samps$FamilyNumber_los), ]
  methylation <- betas[rownames(samps),]
  mydata <- mydata[rownames(samps),]
  
  table(rownames(mydata)==rownames(methylation))

  # Run GEE model
  geemodel <- function(data) {
    gee(
      CpGi ~ tester + smoke + gen_meds + age + sex + NK_Perc + Mono_Perc + Bcell_Perc +
        CD4T_Perc + CD8T_Perc + Array_rownum + Sample_Plate,
      data = data, id = familynumber, family = gaussian,
      corstr = "exchangeable", maxiter = 100, na.action = na.omit
    )
  }
  
  # Testrun for 1 cpg to determine the number of estimates in output (may vary depending on the number of plates on which the individuals were measured for which you have phenotype data)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  mydata$CpGi <- methylation[,1]
  coeff <- summary(geemodel(mydata))$coefficients
  # this the GEE output for our first CpG site:
  coeff
  # these are the results that we will extract:
  coeff[,c(1,4,5)]
  # note that there is no p-value. We will compute the p-value ourselves, like this:
  pval <-  2*pnorm(-abs(coeff[,5]))
  pval
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # Create objects to save the output from gee
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  Nestimates <- nrow(coeff)
  Ncpgs <-ncol(methylation)
  Ncpgs
  pval <- matrix(NA,Ncpgs,Nestimates)
  estimate <-  matrix(NA,Ncpgs,Nestimates)
  RobustSE <-  matrix(NA,Ncpgs,Nestimates)
  RobustZ  <-  matrix(NA,Ncpgs,Nestimates)
  Nestimates
  
  # Run geemodel for all CpGs and save output
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  # Run entire part until the next # at once
  for (i in 1:Ncpgs)
  {
    mydata$CpGi <- methylation[,i]
    coeff <- summary(geemodel(mydata))$coefficients
    for (j in 1: Nestimates)
    {
      pval[i,j] <- coeff[j,4]
      estimate[i,j] <- coeff[j,1]
      RobustSE[i,j] <- coeff[j,4]
      RobustZ[i,j] <- coeff[j,5]
      pval[i,j] <-  2*pnorm(-abs(coeff[j,5]))    
    }
  }
  colnames(estimate)  <- rownames(coeff)
  colnames(RobustSE)  <- rownames(coeff)
  colnames(RobustZ)  <- rownames(coeff)
  colnames(pval)  <- rownames(coeff)
  rownames(estimate) <- colnames(methylation)
  rownames(RobustSE) <- colnames(methylation)
  rownames(RobustZ) <- colnames(methylation)
  rownames(pval) <- colnames(methylation)
  #------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # Make a table with the results
  results <-  data.frame(estimate[,"tester"],RobustSE[,"tester"],RobustZ[,"tester"],pval[,"tester"])
  colnames(results) <- c("estimate","RobustSE","RobustZ","pvalue")
  dim(results)
  results <- merge(results, epic_anno, by=0)
  head(results)     # note that a new column with the q-value has been added
  
  results <- subset(results, CHR !="X" & CHR != "Y")
  
  dim(results)
  
  # Ncpgs <-nrow(results)
  
  bonf <- 0.05 / nrow(results)  # bonferroni corrected alpha (significance threshold)
  signifresults <- results[which(results$pval < bonf),]   # select significant CpG sites after Bonferroni correction
  dim(signifresults)
  
  beta <- matrix(results$estimate)
  ses <- matrix(results$RobustSE)
  colnames(beta) <- "beta"
  colnames(ses) <- "SE"
  bc <- bacon(teststatistics =NULL, effectsizes = beta, standarderrors = ses) 
  bc
  inflation <- inflation(bc)
  
  # Save the results
  avg = mean(mydata$tester)
  med = median(mydata$tester)
  rng = range(mydata$tester)
  age_avg = mean(mydata$age)
  age_rng = range(mydata$age)
  male = length(which(mydata$sex==1))
  female = length(which(mydata$sex==2))
  sig_cpgs = nrow(signifresults)
  output = c(biomarker, nrow(mydata), nrow(results), sig_cpgs, inflation,avg, med, rng, age_avg, age_rng, male, female)
  run_summary = rbind(run_summary, output)
  
  output_filename <- paste0("/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/tp2/tp2_", biomarker, "_ewas.RData")
  sig_output_filename <- paste0("/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/tp2/tp2_", biomarker, "_sig_ewas.RData")
  save(signifresults, estimate, RobustSE, RobustZ, pval, file = sig_output_filename)
  save(results, estimate, RobustSE, RobustZ, pval, file = output_filename)
  
  cat("EWAS results for biomarker:", biomarker, "saved.\n\n")
}

colnames(run_summary) <- c("Biomarker", "Samples", "CpGs", "Sig CpGs", "Mean", "Median", "Min", "Max", "Average Age", "Min Age", "Max_Age", "Males", "Females")
write.csv(run_summary, "/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/tp2/outliers_run_summary_tp2.csv")

# Clean up workspace
rm(list = ls(all = TRUE))

