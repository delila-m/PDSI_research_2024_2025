# Libraries
#####
library(tidyverse)
library(dplyr)
library(ranger)
library(randomForest)
library(caret)
# load in drought functions
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Functions/")
source("drought_functions.R")
#####

# Loading in Previously Cleaned data
#####
# Data Binned to 1st degree
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/CleanedUS_1")
allStatesDf <- readRDS("CleanedUS_1.RDS")

# Load in data binned to 0.5 degree
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
allStatesDf_0.5 <- readRDS("FactorUS_0.5.RDS")

# Previous categorical model for comparison 
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
load("RFAnalysis0.5_factor_updated.Rdata")
#####

# Previous Model results 
#####
# classification metrics for initial model binned to 0.5 degree
sensitivity(test_0.5$USDM_factor, test_0.5$predicted)
# D0        D1        D2        D3        D4        None 
# 0.6473022 0.3637835 0.4383925 0.4192772 0.4340278 0.7061512

confusionMatrix <- confusionMatrix((test_0.5$predicted_factor), (test_0.5$USDM_factor))
summary(confusionMatrix)
# Confusion Matrix (6 x 6) 
# ===
#   D0     D1     D2     D3     D4   None
# D0   324222  29227   8517   1404    161 137351
# D1    58367  58225  15108   3021    313  25020
# D2    11844  17684  33064   5699    650   6480
# D3     2240   4233   8265  12715   1105   1768
# D4      303    702   1214   1447   3125    409
# None 131107  10702   3724    594    101 351402
# ===
#   Overall Statistics (micro average)
# - Accuracy:          0.62
# - Balanced Accuracy: 0.50
# - Sensitivity:       0.62
# - Specificity:       0.92
# - Precision:         0.62
#####

# RF classification using a more in depth gridsearch cv
#####
library(superml)

# split into training and testing
train_gsx <- allStatesDf_0.5 %>% sample_frac(0.80)
test_gsx <- anti_join(allStatesDf_0.5, train_gsx)

# selecting the relevant columns 
gridsearchtrain_x <- train_gsx %>% select(c(-USDM_Avg)) 
gridsearchtest_x <- test_gsx %>% select(c(-USDM_Avg))
# cross validate a new RF model 
rf <- RFTrainer$new(classification = 1, 
                    seed = 3)

gridsearch <- GridSearchCV$new(trainer = rf, 
                               parameters = list(n_estimators = c(50, 100, 200), 
                                                 max_depth = c(2, 3, 6)), 
                               n_folds = 3, 
                               scoring = c('precision'))
# fit model 
gridsearch$fit(gridsearchtrain_x, "USDM_factor")

# look at the best iteration
gridsearch$best_iteration()
# $n_estimators
# [1] 200
# 
# $max_depth
# [1] 2
# 
# $precision_avg
# [1] 0.9854099
# 
# $precision_sd
# [1] 0.001046917

# predict using testing set
gsx_preds <- rf$predict(df = test_gsx)

# look at results 
gsx_test <- as.factor(as.integer(test_gsx$USDM_factor) - 1)
confusionMatrix(gsx_preds, gsx_test)
# Confusion Matrix and Statistics
#####
# Reference
# Prediction      0      1      2      3      4      5
# 0 466758  16222    448     39      0  25601
# 1   8465 140075   5550    157      5     50
# 2     48   3252  68572   2035     26      2
# 3      8     34   1110  28017    413      0
# 4      0      5      8    200   6808      0
# 5  24328    154      7      1      0 473119
# 
# Overall Statistics
# 
# Accuracy : 0.9307          
# 95% CI : (0.9302, 0.9311)
# No Information Rate : 0.3929          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.8965          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: 0 Class: 1 Class: 2 Class: 3 Class: 4 Class: 5
# Sensitivity            0.9343   0.8769  0.90590  0.92013 0.938776   0.9486
# Specificity            0.9452   0.9872  0.99552  0.99874 0.999832   0.9683
# Pos Pred Value         0.9169   0.9078  0.92746  0.94710 0.969662   0.9508
# Neg Pred Value         0.9569   0.9824  0.99405  0.99804 0.999649   0.9669
# Prevalence             0.3929   0.1256  0.05953  0.02395 0.005703   0.3923
# Detection Rate         0.3671   0.1102  0.05393  0.02203 0.005354   0.3721
# Detection Prevalence   0.4004   0.1214  0.05815  0.02327 0.005522   0.3914
# Balanced Accuracy      0.9397   0.9320  0.95071  0.95943 0.969304   0.9584

# save new model 
save(gridsearch, file = "C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/RFAnalysis_0.5_gcv_updated.Rdata")
load("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/RFAnalysis_0.5_gcv_updated.Rdata")
####

# RF model attempting to recreate results from the best iteration of Grid Search 
#####
# using the same training and testing sets from the gridsearch 
rf_recreation <- ranger(USDM_factor ~ PDSI_Avg + bin.x + bin.y,
                        data = train_gsx, 
                        num.trees = 200, 
                        mtry = 2, 
                        classification = TRUE, 
                        verbose = TRUE, 
                        local.importance = TRUE)
rf_predictions <- predict(rf_recreation, test_gsx)

confusionMatrix(rf_predictions$predictions, test_gsx$USDM_factor)

# XGBoost model
#####
# split into training and testing 
train_0.5 <- allStatesDf_0.5 %>% sample_frac(0.80)
test_0.5 <- anti_join(allStatesDf_0.5, train_0.5)

# Get vectors of just the x and y we are training / testing on 
train_x = data.matrix(train_0.5 %>% select(-c(USDM_Avg, bin.x, bin.y)))
train_y = train_0.5$USDM_factor

# Define predictor and response variables in testing set
test_x = data.matrix(test_0.5 %>% select(-c(USDM_Avg, bin.x, bin.y)))
test_y = test_0.5$USDM_factor

## Using caret for categorical prediction
train_ctrl <- trainControl(
  method = "cv", 
  number = 3, 
  verboseIter = TRUE, 
  allowParallel = TRUE
)

# Define the tuning grid for classification
grid_tune <- expand_grid(
  nrounds = c(100, 500, 1000), 
  max_depth = c(2, 4, 6), 
  eta = 0.3, 
  gamma = 0, 
  colsample_bytree = 1, 
  min_child_weight = 1, 
  subsample = 1
)

# Use train() with categorical target
xgb_tune <- train(
  x = train_x, 
  y = train_y
  trControl = train_ctrl, 
  tuneGrid = grid_tune, 
  method = "xgbTree", 
  verbose = TRUE,
  objective = "multi:softprob"  # For multiclass classification
)

# View best parameters
xgb_tune$bestTune

# Make predictions
xgb_pred_factors <- predict(xgb_tune, newdata = test_x)
test_0.5$predicted_factor <- xgb_pred_factors

# save XGBoost object 
save(train_0.5, test_0.5, xgb_tune, file = "C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/XGBAnalysis_0.5_9-29.Rdata")

load("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/XGBAnalysis_0.5.Rdata")
# look at the results 
confusionMatrix <- confusionMatrix(test_0.5$predicted_factor, test_0.5$USDM_factor)
confusionMatrix
#####

# RF model subset by year using annual PDSI data, matching the frequency of the LBDA
#####

# subset the data to the annual summer average
annual_PDSI_0.5 <- summer.average(allStatesDf_0.5)

# create categorical column for USDM measurements
annual_PDSI_0.5 <- usdm.factor(annual_PDSI_0.5)

# Get unique years and sort them
years <- sort(unique(annual_PDSI_0.5$year))

# Calculate how many years make up approximately 20%
num_test_years <- ceiling(length(years) * 0.2)

# Use the most recent consecutive years as test set
test_years <- tail(years, num_test_years)
train_years <- head(years, length(years) - num_test_years)

# Create train and test datasets
train.yearsplit.0.5.factor <- annual_PDSI_0.5[annual_PDSI_0.5$year %in% train_years, ]
test.yearsplit.0.5.factor <- annual_PDSI_0.5[annual_PDSI_0.5$year %in% test_years, ]


# now we can create the new rf model(s)
rf_recreation <- ranger(USDM_Factor ~ PDSI_Avg + bin.x + bin.y,
                        data = train.yearsplit.0.5.factor, 
                        num.trees = 300, 
                        mtry = 2, 
                        classification = TRUE, 
                        verbose = TRUE, 
                        local.importance = TRUE)

rf_predictions <- predict(rf_recreation, test.yearsplit.0.5.factor)
test.yearsplit.0.5.factor$predictions <- rf_predictions$predictions
confusionMatrix(rf_predictions$predictions, test.yearsplit.0.5.factor$USDM_Factor)
# Confusion Matrix and Statistics
#####
#             Reference
# Prediction   D0 None   D1   D2   D3   D4
# D0         3890 1255  831  235   46    6
# None       4128 3983  405   61    8    0
# D1         458   80  525  307  132   12
# D2         81   17  221  284  185   50
# D3         8    2   53  114  169   21
# D4         1    0    5   17   47   23
# 
# Overall Statistics
# 
# Accuracy : 0.5025          
# 95% CI : (0.4951, 0.5099)
# No Information Rate : 0.4851          
# P-Value [Acc > NIR] : 1.839e-06       
# 
# Kappa : 0.2549          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: D0 Class: None Class: D1 Class: D2 Class: D3 Class: D4
# Sensitivity             0.4541      0.7463   0.25735   0.27898   0.28790  0.205357
# Specificity             0.7391      0.6266   0.93668   0.96671   0.98840  0.996011
# Pos Pred Value          0.6211      0.4639   0.34676   0.33890   0.46049  0.247312
# Neg Pred Value          0.5897      0.8508   0.90617   0.95637   0.97583  0.994934
# Prevalence              0.4851      0.3022   0.11552   0.05764   0.03324  0.006342
# Detection Rate          0.2203      0.2255   0.02973   0.01608   0.00957  0.001302
# Detection Prevalence    0.3546      0.4861   0.08573   0.04745   0.02078  0.005266
# Balanced Accuracy       0.5966      0.6864   0.59702   0.62284   0.63815  0.600684
#####

# save model for future use 
save(rf_recreation, annual_PDSI_0.5, train.yearsplit.0.5.factor, test.yearsplit.0.5.factor, rf_predictions, 
     file = "C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/rf.annual.yearsplit.0.5.Rdata")
load("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/rf.annual.yearsplit.0.5.Rdata")

# plot results for 2024 
preds.2024 <- test.yearsplit.0.5.factor %>% filter(year == 2024)
plotdata <- pred.v.actual.plot.factor(preds.2024, "USDM_Factor", "predictions", 
                          save = FALSE, "testplot", year = 2024)







