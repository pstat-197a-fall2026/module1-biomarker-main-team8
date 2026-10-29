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

bench <- readRDS('results/benchmark.rds')

bench$s1      # in-class t-test top 10
bench$s2      # in-class random forest top 10
bench$sstar   # in-class panel (intersection)
bench$metrics # in-class sens, spec, acc, AUROC



proteins_s1_train
bench$s1
intersect(bench$s1, proteins_s1_train)   # kept
setdiff(bench$s1, proteins_s1_train)     # dropped
setdiff(proteins_s1_train, bench$s1)     # new

sum(ttests_train$p.adj < 0.05)

# Q3a t-tests (training only) vs in-class top 10:
#   kept 8/10: DERM, RELT, Calcineurin, C1QR1, MRC2, IgD, PTN, FSTL1
#   dropped: CXCL16, Cadherin-5 | new: TSP4, RET


#   only DERM still has p.adj < 0.05 (less power with 123 children)
#   -> t-test selection fairly stable, but much weaker significance


