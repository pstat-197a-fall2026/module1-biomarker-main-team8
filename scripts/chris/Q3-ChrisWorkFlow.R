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

# checking the split worked properly
biomarker_split
table(train$group)
table(test$group)



## STEP 2: MULTIPLE TESTING (training data only)
#################################################

# function to compute tests (same as class)
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_train <- train %>%                 # <- train, not biomarker_clean
  select(-ados) %>%
  pivot_longer(-group,
               names_to = 'protein',
               values_to = 'level') %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# top 10 proteins by adjusted p-value (training data only)
proteins_s1_train <- ttests_train %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

# compare with the in-class top 10
proteins_s1_train
sum(ttests_train$p.adj < 0.05)