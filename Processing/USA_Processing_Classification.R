# Libraries
library(tidyverse)
library(dplyr)
library(ranger)

# Data Binned to 1st degree
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/CleanedUS_1")
allStatesDf <- readRDS("CleanedUS_1.RDS")


# Load in data binned to 0.5 degree
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
allStatesDf_0.5 <- readRDS("FactorUS_0.5.RDS")

setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
load("RFAnalysis0.5_factor_updated.Rdata")

# classification metrics for model binned to 0.5 degree

sensitivity(test_0.5$USDM_factor, test_0.5$predicted)
# D0        D1        D2        D3        D4        None 
# 0.6473022 0.3637835 0.4383925 0.4192772 0.4340278 0.7061512

confusionMatrix <- cmatrix(test_0.5$USDM_factor, test_0.5$predicted)
summary(confusionMatrix)
# Confusion Matrix (6 x 6) 
# ===============================================================
#   D0     D1     D2     D3     D4   None
# D0   324222  29227   8517   1404    161 137351
# D1    58367  58225  15108   3021    313  25020
# D2    11844  17684  33064   5699    650   6480
# D3     2240   4233   8265  12715   1105   1768
# D4      303    702   1214   1447   3125    409
# None 131107  10702   3724    594    101 351402
# ===============================================================
#   Overall Statistics (micro average)
# - Accuracy:          0.62
# - Balanced Accuracy: 0.50
# - Sensitivity:       0.62
# - Specificity:       0.92
# - Precision:         0.62


# classification random forest model 
# split into training and testing 
train <- grouped %>% sample_frac(0.80)
test <- anti_join(grouped , train)

# trying rf using a more indepth gridsearch cv
install.packages("superml")
library(superml)

rf <- RFTrainer$new()
gridsearch <- GridSearchCV$new(trainer = rf, 
                               parameters = list(








