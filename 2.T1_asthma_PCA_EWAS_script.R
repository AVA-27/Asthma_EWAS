# Asthma PCA and EWAS Script

rm(list = ls())


# Load packages
library(rospca)
library(DNAmArray)
library(gee)
library(ggplot2)
library(dplyr)
library(pheatmap)
# library(FactoMineR)
# library(factoextra)


# Load methylation data and sample/individual information
load("/data/avavanasselt/data/betas_final.Rdata")
load("/data/avavanasselt/sampsheets/sampsheet_OLRM_2023.RData")

# Perform same selection of samples as done in script 1
# Select TP1 samples
tp1blood <- sampsheet[which(sampsheet$TP1_sample == 1), ]
tp1blood <- tp1blood[order(tp1blood$FamilyNumber_los), ]
dim(tp1blood)
# [1] 341 282
# remove 3 re-run samples
dupes <- tp1blood[which(duplicated(tp1blood$FISNumber)), ]
dim(dupes)
# [1]   6 282

tp1blood <- tp1blood[-which(duplicated(tp1blood$FISNumber)), ]
dim(tp1blood)
# [1] 335 282

table(tp1blood$smoke, tp1blood$smoke_quit)
# -9   0   1
# -9   1   0   0
# 0    9 177  61
# 1    2  85   0
# remove people with both missing status
dim(tp1blood)
# [1] 335 282
tp1blood <- subset(tp1blood, smoke != -9 & smoke_quit != -9)
dim(tp1blood)
# [1] 323 282

table(tp1blood$smoke, tp1blood$smoke_quit)
# 0   1
# 0 177  61
# 1  85   0


# create new column with updated smoking status (0 = never smoke, 1 = previous smoker, 2 = current smoker)
tp1blood <- tp1blood %>%
  mutate(smoke_status = if_else(smoke == 0 & smoke_quit == 1, 1, if_else(smoke == 1, 2, 0)))


table(tp1blood$smoke_status)
# 0   1   2 
# 177  61  85 
# Subset betas based on TP1 samples
betas <- betas[, rownames(tp1blood)]

# Transpose betas
betas <- t(betas)

table(rownames(tp1blood)==rownames(betas))
# TRUE 
# 323 


# load imputed data from script 1
load("/data/avavanasselt/results/asthma/pca/mice_impute/pheno_imputed_mice.RData")
pheno_data <- as.data.frame(pheno_imputed)
dim(pheno_data)
#[1] 319  35 

#convert data type for scaling and PCA
pheno_data <- pheno_data %>% mutate_if(is.factor, as.numeric)

# look at data structure
means <- data.frame(sapply(pheno_data, mean))
variance <- data.frame(sapply(pheno_data, var))
sd <- data.frame(sapply(pheno_data, sd))

stats <- merge(means, variance, by=0)
rownames(stats) <- stats$Row.names
stats <- merge(stats, sd, by=0)
stats <- subset(stats, select= -c(1))
colnames(stats) <- c("Variables", "Mean", "Variance", "Standard Deviation")
write.csv(stats, "/data/avavanasselt/results/asthma/pca/variable_descriptors.csv")


### SCALING ###
# z-score transformation - treat all variables as continuous
scaled_pheno <- as.data.frame(scale(pheno_data))



### PCA ###
# Set up PCA analysis
pca <- robpca(scaled_pheno, k = 10, kmax = 10, alpha = 0.75, h = NULL, mcd = FALSE,
        ndir = "all", skew = TRUE)


pca_values <- pca$loadings


# $eigenvalues
# PC1       PC2       PC3       PC4       PC5       PC6       PC7       PC8 
# 3.6020754 2.3002297 1.5325243 1.2480955 0.7212765 0.4841179 0.3234200 0.1251906 
# PC9      PC10 
# 0.1162162 0.0542785 

# sd value in pca 
# 13.69292

# sum(pca$eigenvalues)
# [1] 12.71254

# sum(pca$eigenvalues)/13.69292
# [1] 0.9284024


### Save all data from PCA ###
pdf("/data/avavanasselt/results/asthma/pca/scaled_pca_phenotype_heatmap.pdf")
pheatmap(pca_values, cluster_rows=TRUE, cluster_cols = FALSE, color=colorRampPalette(c("coral", "azure", "aquamarine"))(50))
dev.off()

save(pca_values, file = "/data/avavanasselt/results/asthma/pca/pca_phenotype_values.RData")
pca_scores <- pca$scores
save(pca_scores, file = "/data/avavanasselt/results/asthma/pca/pca_sample_values.RData")
sum(pca$eigenvalues)
eigens <- pca$eigenvalues
eigens
save(eigens, file = "/data/avavanasselt/results/asthma/pca/pca_eigenvalues.RData")

pca$scores
# object for each PCA for each sample/individual


# # remove outliers
library(ewaff)

pca_scores <- t(pca_scores)
pca_scores <- ewaff.handle.outliers(pca_scores, method="iqr", iqr.limit=3)[[1]]
pca_scores <- t(pca_scores)

# change pcs 8-10 to Y/N valiable instead
pca_scores <- data.frame(pca_scores)

pca_scores = pca_scores %>% 
  mutate(PC8_out = if_else(is.na(PC8), 1, 0))

pca_scores = pca_scores %>% 
  mutate(PC9_out = if_else(is.na(PC9), 1, 0))

pca_scores = pca_scores %>% 
  mutate(PC10_out = if_else(is.na(PC10), 1, 0))

save(pca_scores, file = "/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/pca_scores_final.RData")



####  Run EWAS  #######
npcas <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8_out", "PC9_out", "PC10_out")

#Order TP1 data 
tp1blood <- tp1blood[rownames(pheno_data),]
dim(tp1blood)
betas <- betas[rownames(tp1blood),]

table(rownames(tp1blood)==rownames(pheno_data))
table(rownames(tp1blood)==rownames(betas))
table(rownames(pca_scores)==rownames(tp1blood))


# Create a data frame for analysis variables
data <- data.frame(
  familynumber = as.numeric(tp1blood$FamilyNumber_los),
  CpGi = rep(NA, nrow(betas)),
  peak_flow = as.numeric(tp1blood$peak_flow),
  sex = as.numeric(tp1blood$sex_los),
  age = as.numeric(tp1blood$age),
  NK_Perc = as.numeric(tp1blood$NK),
  Mono_Perc = as.numeric(tp1blood$Mono),
  Bcell_Perc = as.numeric(tp1blood$Bcell),
  CD4T_Perc = as.numeric(tp1blood$CD4T),
  CD8T_Perc = as.numeric(tp1blood$CD8T),
  Array_rownum = as.numeric(tp1blood$Array_rownum),
  Sample_Plate = as.factor(tp1blood$Box_plate),
  height = as.numeric(tp1blood$height),
  weight = as.numeric(tp1blood$weight),
  smoke = as.numeric(tp1blood$smoke_status),
  med_lungs = as.numeric(tp1blood$med_lungs),
  betamimetica1 = as.numeric(tp1blood$betamimetica1),
  inh_corticosteroid1 = as.numeric(tp1blood$inh_corticosteroid1)
)
rownames(data) <- rownames(tp1blood)


# Make dataframe to record run summary data for each EWAS
run_summary = data.frame()

# Iterate through biomarkers and run EWAS
for (biomarker in npcas) {
  cat("Running EWAS for biomarker:", biomarker, "\n")
  data$tester <- pca_scores[,biomarker]
  mydata <- data[complete.cases(data[, "tester"]), ]
  samps <- tp1blood[rownames(mydata),]
  samps <- samps[order(samps$FamilyNumber_los), ]
  methylation <- betas[rownames(samps),]
  mydata <- mydata[rownames(samps),]
  
  table(rownames(mydata)==rownames(methylation))
  # Run GEE model
  geemodel <- function(data) {
    gee(
      CpGi ~ tester + smoke + med_lungs + betamimetica1 + inh_corticosteroid1 + age + sex + NK_Perc + Mono_Perc + Bcell_Perc +
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
  
  bonf <- 0.05 / nrow(results)  # bonferroni corrected alpha (significance threshold)
  signifresults <- results[which(results$pval < bonf),]   # select significant CpG sites after Bonferroni correction
  dim(signifresults)
  
  # Save the results
  avg = mean(mydata$tester)
  med = median(mydata$tester)
  rng = range(mydata$tester)
  age_avg = mean(mydata$age)
  age_rng = range(mydata$age)
  male = length(which(mydata$sex==1))
  female = length(which(mydata$sex==2))
  sig_cpgs = nrow(signifresults)
  output = c(biomarker, nrow(mydata), Ncpgs, sig_cpgs, avg, med, rng, age_avg, age_rng, male, female)
  run_summary = rbind(run_summary, output)
  
  output_filename <- paste0("/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/", biomarker, "_ewas.RData")
  sig_output_filename <- paste0("/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/", biomarker, "_sig_ewas.RData")
  save(signifresults, estimate, RobustSE, RobustZ, pval, file = sig_output_filename)
  save(results, estimate, RobustSE, RobustZ, pval, file = output_filename)
  
  cat("EWAS results for biomarker:", biomarker, "saved.\n\n")
}

# Note - we removed sex chromosome data in a later script and re-calculated significane threshold at that time

colnames(run_summary) <- c("Biomarker", "Samples", "CpGs", "Sig CpGs", "Mean", "Median", "Min", "Max", "Average Age", "Min Age", "Max_Age", "Males", "Females")
write.csv(run_summary, "/data/avavanasselt/results/asthma/pca/mice_impute/olrm_medcovs/run_summary_PCA.csv")

# Clean up workspace
rm(list = ls(all = TRUE))

