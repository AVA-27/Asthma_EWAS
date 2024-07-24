# SCRIPT MADE FOR OBSERVING CLINICAL MARKER DATA AND IMPUTING MISSING VALUES #

rm(list = ls())

# Load packages 
library(rospca)
library(DNAmArray)
library(gee)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(SparseArray)
library(cellWise)
library(mice)
library(cowplot)


# Load methylation data and sample/individual information
load("/data/avavanasselt/data/betas_final.Rdata")
load("/data/avavanasselt/sampsheets/sampsheet_OLRM_2023.RData")

# Select TP1 samples
tp1blood <- sampsheet[which(sampsheet$TP1_sample == 1), ]
tp1blood <- tp1blood[order(tp1blood$FamilyNumber_los), ]
dim(tp1blood)
# [1] 341 282

# remove 3 re-run samples
dupes <- tp1blood[which(duplicated(tp1blood$FISNumber)), ]
dim(dupes)
tp1blood <- tp1blood[-which(duplicated(tp1blood$FISNumber)), ]
dim(tp1blood)
# [1] 335 282


# Check and update smoking status
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


# Gather clinical marker data that we wish to measure
pheno_data <- select(tp1blood, bronch_hyperreactivity2, atopy, reversibility, 
                     peak_flow, vital_capacity, forc_exp_vol, tiffenau_index, pc20m, pc10m, eosin_granulocyt, immunoglobulin,
                     dustmites_IgE, cats_IgE, grass_IgE, coughing, wheezing, dysp, dyspn, bronch_hyperreactivity1, 
                     med_lungs, betamimetica1, inh_corticosteroid1,skin_reaction_degr, grasses, trees1,trees2, 
                     weeds, mites1, mites2, cat, dog, horse, hairs, feathers, asp.fum, clad.h, alt.alte, cand.al)

dim(pheno_data)
# [1] 323  38

############## count number of -9s per row (-9 can indicate missing)

#first convert N/As and -8's to -9 to match rest of data set
colSums(is.na(pheno_data))
# 5 instances in peak flow

# convert NAs of -9
pheno_data[is.na(pheno_data)] <- -9


test <- apply(pheno_data, 2, function(x) length(which(x == -8)))
# 18 instances in PC10m

# convert to -9
pheno_data <- replace(pheno_data, pheno_data == -8, -9)

ind_counts <- apply(pheno_data, 1, function(x) length(which(x==-9)))

table(ind_counts)
#  0   1   2   3   4   5   6  17  20  22 
# 186  68  27  14  20   1   3   2   1   1 


biomarker_counts <- apply(pheno_data, 2, function(x) length(which(x==-9)))

table(biomarker_counts)
#  0  1  2  4  7 11 14 20 26 29 32 48 59 
# 10  2  1 16  1  1  1  1  1  2  1  1  1 


# remove samples with more than 6 missing values
pheno_data <- subset(pheno_data, ind_counts < 7)
dim(pheno_data)
# [1] 319  38


ind_counts_2 <- apply(pheno_data, 1, function(x) length(which(x==-9)))

table(ind_counts_2)


biomarker_counts_2 <- apply(pheno_data, 2, function(x) length(which(x==-9)))

table(biomarker_counts_2)

# convert -9 tp NA values

pheno_data[pheno_data == -9] <- NA

colSums(is.na(pheno_data))

table(colSums(is.na(pheno_data)))

# Remove medication data from clinical marker data (will be used as covariates later)
pheno_data$med_lungs <- NULL
pheno_data$betamimetica1 <- NULL
pheno_data$inh_corticosteroid1 <- NULL

dim(pheno_data)
# [1] 319  35

# Prep data for imputation by modifying data type where necessary
str(pheno_data)
# convert "int" to "factor" where needed
pheno_data$bronch_hyperreactivity2 = as.factor(pheno_data$bronch_hyperreactivity2)
pheno_data$atopy = as.factor(pheno_data$atopy)
pheno_data$reversibility = as.factor(pheno_data$reversibility)
# leave - $ peak_flow              : num  10.2 7.3 8.91 9.1 9.51 ...
# leave - $ vital_capacity         : num  4.87 4.46 4.23 4.15 5.68 4.16 4.35 3.75 6.2 4.83 ...
# leave - $ forc_exp_vol           : num  4.12 3.93 3.25 3.66 4.64 3 3.32 3.52 4.96 4.38 ...
# leave - $ tiffenau_index         : num  0.85 0.88 0.77 0.88 0.82 0.72 0.76 0.93 0.8 0.91 ...
#  included - leave - $ pc20m                  : num  160 160 160 160 80 140 5 160 160 160 ...
#  included - leave - $ pc10m                  : num  160 -8 7.5 -8 50 50 3.12 160 120 35 ...
#leave - $ eosin_granulocyt       : num  0.11 0.06 0.07 0.08 0.16 0.18 0.14 0.14 0.07 0.05 ...
pheno_data$immunoglobulin = as.numeric(pheno_data$immunoglobulin)
# leave - $ dustmites_IgE          : num  0 0 0 0 0 0 0 0 0 0 ...
# leave - $ cats_IgE               : num  0 0 68.2 0 0 0 0 0 0 0 ...
# leave - $ grass_IgE              : num  0 0 4.2 0 96.1 0 0 0 0 0 ...
pheno_data$coughing = as.factor(pheno_data$coughing)
pheno_data$wheezing = as.factor(pheno_data$wheezing)
pheno_data$dysp  = as.factor(pheno_data$dysp)
pheno_data$dyspn = as.factor(pheno_data$dyspn)
pheno_data$bronch_hyperreactivity1 = as.factor(pheno_data$bronch_hyperreactivity1)
# not used - pheno_data$med_lungs = as.factor(pheno_data$med_lungs)
# not used - pheno_data$betamimetica1 = as.factor(pheno_data$betamimetica1)
# not used - pheno_data$inh_corticosteroid1 = as.factor(pheno_data$inh_corticosteroid1)
pheno_data$skin_reaction_degr = as.factor(pheno_data$skin_reaction_degr)
pheno_data$grasses = as.factor(pheno_data$grasses)
pheno_data$trees1 = as.factor(pheno_data$trees1)
pheno_data$trees2 = as.factor(pheno_data$trees2)
pheno_data$weeds = as.factor(pheno_data$weeds)
pheno_data$mites1 = as.factor(pheno_data$mites1)
pheno_data$mites2 = as.factor(pheno_data$mites2)
pheno_data$cat = as.factor(pheno_data$cat)
pheno_data$dog = as.factor(pheno_data$dog)
pheno_data$horse = as.factor(pheno_data$horse)
pheno_data$hairs = as.factor(pheno_data$hairs)
pheno_data$feathers = as.factor(pheno_data$feathers)
pheno_data$asp.fum = as.factor(pheno_data$asp.fum)
pheno_data$clad.h = as.factor(pheno_data$clad.h)
pheno_data$alt.alte = as.factor(pheno_data$alt.alte)
pheno_data$cand.al = as.factor(pheno_data$cand.al)


# Perform imputation
pheno_impute <- mice(pheno_data, seed = 27)
pheno_imputed <- data.frame(complete(mice(pheno_data, seed = 27)))

# Plot a grid of histograms and see what it looks like
h1 <- ggplot(pheno_data, aes(x = bronch_hyperreactivity2)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity", stat = "count") +
  ggtitle("BH2 Original distribution") +
  theme_classic()
h2 <- ggplot(pheno_imputed, aes(x = bronch_hyperreactivity2)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity", stat = "count") +
  ggtitle("BH2 Imputed distribution") +
  theme_classic()
h3 <- ggplot(pheno_data, aes(x = reversibility)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity", stat = "count") +
  ggtitle("Reversibility Original distribution") +
  theme_classic()
h4 <- ggplot(pheno_imputed, aes(x = reversibility)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity", stat = "count") +
  ggtitle("Reversibility Imputed distribution") +
  theme_classic()
h5 <- ggplot(pheno_data, aes(x = peak_flow)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Peak Flow Original distribution") +
  theme_classic()
h6 <- ggplot(pheno_imputed, aes(x = peak_flow)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Peak Flow Imputed distribution") +
  theme_classic()
h7 <- ggplot(pheno_data, aes(x = forc_exp_vol)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("FEV Original distribution") +
  theme_classic()
h8 <- ggplot(pheno_imputed, aes(x = forc_exp_vol)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("FEV Imputed distribution") +
  theme_classic()
h9 <- ggplot(pheno_data, aes(x = tiffenau_index)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("TI Original distribution") +
  theme_classic()
h10 <- ggplot(pheno_imputed, aes(x = tiffenau_index)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("TI Imputed distribution") +
  theme_classic()
h11 <- ggplot(pheno_data, aes(x = eosin_granulocyt)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Eosinophil Original distribution") +
  theme_classic()
h12 <- ggplot(pheno_imputed, aes(x = eosin_granulocyt)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Eosinophil Imputed distribution") +
  theme_classic()
h13 <- ggplot(pheno_data, aes(x = immunoglobulin)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Immunoglobulin Original distribution") +
  theme_classic()
h14 <- ggplot(pheno_imputed, aes(x = immunoglobulin)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Immunoglobulin Imputed distribution") +
  theme_classic()
h15 <- ggplot(pheno_data, aes(x = dustmites_IgE)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Dustmites IgE Original distribution") +
  theme_classic()
h16 <- ggplot(pheno_imputed, aes(x = dustmites_IgE)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Dustmites IgE Imputed distribution") +
  theme_classic()
h17 <- ggplot(pheno_data, aes(x = cats_IgE)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Cats IgE Original distribution") +
  theme_classic()
h18 <- ggplot(pheno_imputed, aes(x = cats_IgE)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Cats IgE Imputed distribution") +
  theme_classic()
h19 <- ggplot(pheno_data, aes(x = grass_IgE)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Grass IgE Original distribution") +
  theme_classic()
h20 <- ggplot(pheno_imputed, aes(x = grass_IgE)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Grass IgE Imputed distribution") +
  theme_classic()
h21 <- ggplot(pheno_data, aes(x = pc10m)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Pc10m Original distribution") +
  theme_classic()
h22 <- ggplot(pheno_imputed, aes(x = pc10m)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Pc10m Imputed distribution") +
  theme_classic()
h23 <- ggplot(pheno_data, aes(x = pc20m)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Pc20m Original distribution") +
  theme_classic()
h24 <- ggplot(pheno_imputed, aes(x = pc20m)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Pc20m Imputed distribution") +
  theme_classic()

png(file = paste0("/data/avavanasselt/results/asthma/pca/mice_impute/imputed_distributions.png"), height = 80, width = 15, units="in", res=300)
plot_grid(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, nrow = 12, ncol = 2)
dev.off()

# save data for future analysis
save(pheno_imputed, file = "/data/avavanasselt/results/asthma/pca/mice_impute/pheno_imputed_mice.RData")

save(pheno_data, file = "/data/avavanasselt/results/asthma/pca/pheno_data.RData")


