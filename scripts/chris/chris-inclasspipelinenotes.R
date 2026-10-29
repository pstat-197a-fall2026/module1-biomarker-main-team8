# ---------------------------------------------------------------------
# Data:    154 children (76 ASD, 78 TD) x 1,317 proteins;
#          values = log10 level, z-scored across all 154, clipped at +/-3
# Steps:   t-tests (top 10) + random forest (top 10) -> intersection
#          -> logistic regression on 80% -> scored on 20% (31 children)
# Panel:   DERM, RELT, MRC2, IgD, Cadherin-5
# Results: sens 0.812, spec 0.733, acc 0.774, AUROC 0.883
# Caveat:  proteins chosen using all 154 children (incl. test set) and
#          one small split, so results are likely optimistic (see Q3a)
# =====================================================================

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group,
               names_to = 'protein',
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>%
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)
# Welch t-test per protein: H0 mean ASD = mean TD, two-sided
# p.adj = Benjamini-Yekutieli (controls false discovery rate)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

# Top 10 all lower in ASD (0.75-0.89 SD); DERM largest gap (-0.885)
# All 10 significant after BY correction (p.adj 0.00008-0.0018)
# Note: p.adj order differs from raw p order (FSTL1 3rd by p, 9th by p.adj),
#       but the same 10 proteins are selected at k = 10

sum(ttests_out$p.adj < 0.05)
# 113 of 1,317 proteins pass at p.adj < 0.05

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))
# 154 x 1,317 table: each entry = SDs from that protein's overall average

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors,
                       y = response,
                       ntree = 1000,
                       importance = T)
# 1,000 trees, 36 proteins tried at each split (default sqrt(1,317))

# check errors
rf_out$confusion
# OOB error 31.2%: ASD 35.5% misclassified (27/76), TD 26.9% (21/78)

# compute importance scores
proteins_s2 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)
# Top 10 by Gini importance: DERM, IgD, MAPK14, RELT, eIF-4H,
# TSP4, MRC2, ALCAM, FSTL3, Cadherin-5 (changes with seed/mtry)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)
# proteins_s1 from multiple testing, proteins_s2 from random forest;
# keep proteins in both lists
# Panel: DERM, RELT, MRC2, IgD, Cadherin-5

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)
# class: TRUE = ASD, FALSE = TD

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)
# 123 training / 31 test children (split made AFTER protein selection)

# fit logistic regression model to training set
fit <- glm(class ~ .,
           data = training(biomarker_split),
           family = 'binomial')
# log-odds of ASD = b0 + b1*DERM + b2*RELT + b3*MRC2 + b4*IgD + b5*Cadherin-5
# Deviance 170.4 -> 131 on 5 df (panel beats intercept-only model)

# evaluate errors on test set
class_metrics <- metric_set(sensitivity,
                            specificity,
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred,
                event_level = 'second')
# Factors built in mutate() first (newer yardstick rejects factor() inside)
# Benchmark: sens 0.812, spec 0.733, acc 0.774, AUROC 0.883 (31 children)

summary(fit)        # Estimate, Std. Error, z value, Pr(>|z|) per protein
broom::tidy(fit)    # the same as a tibble (tidymodels loads broom)
# Significant: DERM (p = 0.020), IgD (p = 0.004)
# RELT, MRC2, Cadherin-5 not significant (overlap with DERM and IgD)

# odds ratios: e^beta = multiplicative change in odds of ASD per +1 SD of each protein
exp(coef(fit))
# DERM OR 0.47: +1 SD -> ~53% lower odds of ASD
# IgD  OR 0.52: +1 SD -> ~48% lower odds of ASD

# odds ratios with 95% CIs; an interval containing 1 means not significant
broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
















# benchmark row for Q3 comparisons
benchmark <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred, event_level = 'second') %>%
  mutate(variant = 'in-class benchmark', n_proteins = length(proteins_sstar))

saveRDS(benchmark, 'results/benchmark.rds')

saveRDS(list(metrics = benchmark,
             s1 = proteins_s1, s2 = proteins_s2, sstar = proteins_sstar),
        'results/benchmark.rds')


