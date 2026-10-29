# ---- Q3a: selection on a training partition ----
# Split 80/20 FIRST, before any selection
# Run t-tests, random forest and intersection on training data only
# Fit logistic regression on training data; score once on the test data
# Goal: see how much the benchmark was inflated by selecting
#       proteins with all 154 children (incl. the test set)

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')



## STEP 1: SPLIT FIRST
#######################
set.seed(101422)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8)

train <- training(biomarker_split)   # 123 children: used for everything
test  <- testing(biomarker_split)    # 31 children: untouched until the end

# check
biomarker_split
table(train$group)
table(test$group)