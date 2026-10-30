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



## STEP 3: RANDOM FOREST (training data only)
###############################################

predictors_train <- train %>% select(-c(group, ados))
response_train   <- train %>% pull(group) %>% factor()

set.seed(101422)
rf_train <- randomForest(x = predictors_train,
                         y = response_train,
                         ntree = 1000,
                         importance = T)

rf_train$confusion   # out-of-bag errors on training children

proteins_s2_train <- rf_train$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_train$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

# compare with in-class forest top 10
proteins_s2_train
intersect(bench$s2, proteins_s2_train)   # kept
setdiff(bench$s2, proteins_s2_train)     # dropped
setdiff(proteins_s2_train, bench$s2)     # new


# Q3a forest (training only) vs in-class forest top 10:
#   OOB error 32.5% (ASD 43.3%, TD 22.2%) vs 31.2% in class
#   kept 6/10: DERM, IgD, MAPK14, eIF-4H, TSP4, ALCAM
#   dropped: RELT, MRC2, FSTL3, Cadherin-5 | new: CSK, FSTL1, CK-MB, 14-3-3 zeta/delta
#   -> forest list less stable than t-test list; 3 of 5 class panel proteins gone



## STEP 4: PANEL + LOGISTIC REGRESSION (training data only)
############################################################

proteins_sstar_train <- intersect(proteins_s1_train, proteins_s2_train)
proteins_sstar_train
intersect(bench$sstar, proteins_sstar_train)   # kept from class panel

train_sstar <- train %>%
  select(group, any_of(proteins_sstar_train)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

fit_train <- glm(class ~ ., data = train_sstar, family = 'binomial')
broom::tidy(fit_train)


# Q3a panel (training only): DERM, IgD, TSP4, FSTL1 (4 vs 5 in class)
#   only DERM, IgD shared with class panel; both still the only significant terms
#   DERM -0.716 (p = 0.028), IgD -0.645 (p = 0.004), similar to class
#   AIC 140.8 vs 143 in class: similar fit with one fewer protein
#   same seed + same 154 rows -> same 123/31 split as the benchmark

## STEP 5: EVALUATE ON TEST SET (first use of test data)
#########################################################

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

q3a <- test %>%
  mutate(class = (group == 'ASD')) %>%
  add_predictions(fit_train, type = 'response') %>%
  mutate(truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred, event_level = 'second') %>%
  mutate(variant = 'Q3a: split first', n_proteins = length(proteins_sstar_train))

# compare with benchmark
bind_rows(bench$metrics, q3a) %>%
  select(variant, n_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

# Q3a result vs benchmark (same 31 test children):
#   sens 0.750 vs 0.812 | spec 0.800 vs 0.733 | acc 0.774 vs 0.774 | AUROC 0.871 vs 0.883
#   panel changed (3 of 5 proteins swapped) but accuracy nearly unchanged:
#   DERM and IgD drive the model and survive either way
#   -> leakage inflated the benchmark only slightly on this split (AUROC -0.012)
#   caveat: one 31-child split; repeat over many splits to confirm

# sensitivity: share of ASD children correctly identified (ASD caught / all ASD)
# specificity: share of TD children correctly identified (TD cleared / all TD)
# accuracy:    share of all children classified correctly (correct / all)



## ===================== Q3b: LARGER k =====================
# Choose more than 10 top proteins per method (k = 15, 20, 30, 50)
# Reuses the Q3a split-first setup: train/test, ttests_train, rf_train


## STEP 1: PANEL FOR EACH k
############################

# full rankings from the training data (best first)
rank_ttest <- ttests_train %>% arrange(p.adj) %>% pull(protein)

rank_rf <- rf_train$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_train$importance)) %>%
  arrange(desc(MeanDecreaseGini)) %>%
  pull(protein)

# intersection panel for each k
k_values <- c(10, 15, 20, 30, 50)

q3b_panels <- tibble(k = k_values) %>%
  mutate(panel      = map(k, ~ intersect(head(rank_ttest, .x), head(rank_rf, .x))),
         n_proteins = map_int(panel, length))

q3b_panels               # panel size for each k
q3b_panels$panel[[1]]    # k = 10: should match Q3a (DERM, IgD, TSP4, FSTL1)
q3b_panels$panel[[2]]    # k = 15
q3b_panels$panel[[3]]    # k = 20
q3b_panels$panel[[4]]    # k = 30
q3b_panels$panel[[5]]    # k = 50


# Q3b panels: k = 10/15/20/30/50 -> 4/5/6/12/29 proteins
#   k = 10 matches Q3a (DERM, IgD, TSP4, FSTL1); panels are nested
#   slow growth to k = 20, jumps at k = 30 (RELT, MRC2, Cadherin-5 return)
#   k = 50: 29 proteins on 123 children -> overfitting risk

## STEP 2: FIT + SCORE EACH PANEL
###################################

fit_score <- function(panel){
  train_k <- train %>%
    select(group, all_of(panel)) %>%
    mutate(class = (group == 'ASD')) %>%
    select(-group)
  
  fit_k <- glm(class ~ ., data = train_k, family = 'binomial')
  
  test %>%
    mutate(class = (group == 'ASD')) %>%
    add_predictions(fit_k, type = 'response') %>%
    mutate(truth = factor(class, levels = c(FALSE, TRUE)),
           est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
    class_metrics(truth = truth, estimate = est, pred, event_level = 'second')
}

q3b <- q3b_panels %>%
  mutate(metrics = map(panel, fit_score)) %>%
  select(k, n_proteins, metrics) %>%
  unnest(metrics) %>%
  mutate(variant = paste0('Q3b: k = ', k))

# comparison table
bind_rows(bench$metrics, q3a, q3b) %>%
  select(variant, n_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# Q3b results (same 31 test children):
#   k = 10/15/20: 4-6 proteins, acc 0.774-0.839, AUROC 0.858-0.871 (best range)
#   k = 15 best accuracy (0.839) but same AUROC as k = 10 -> 2-child gain, not clear win
#   k = 30/50: 12-29 proteins, acc 0.613-0.710, AUROC 0.758-0.804 (overfitting)
#   -> more proteins does not help; small panels around DERM + IgD do best


q3b %>%
  filter(.metric %in% c('accuracy', 'roc_auc')) %>%
  ggplot(aes(k, .estimate, color = .metric)) +
  geom_line() + geom_point() +
  labs(x = 'k (top proteins per method)', y = 'Test estimate')

# Plot: test accuracy and AUROC vs k (top proteins per method)
#   AUROC flat at 0.871 for k = 10-15, slips at k = 20, drops sharply at k = 30
#   accuracy peaks at k = 15 (0.839), falls to 0.613 at k = 30
#   partial recovery at k = 50 is likely single-split noise, not a real gain
#   -> small panels (k <= 20) do best; large panels overfit

## ===================== Q3c: FUZZY INTERSECTION =====================
# Combine top proteins across methods with a looser rule than intersect()
# Reuses Q3a/Q3b objects: rank_ttest, rank_rf, fit_score, train/test

## FUZZY PANELS
#################################
# Rule A: in the top k of one list AND the top 2k of the other

fuzzy_panel <- function(k){
  union(intersect(head(rank_ttest, k),     head(rank_rf, 2 * k)),
        intersect(head(rank_ttest, 2 * k), head(rank_rf, k)))
}

q3c_panels <- tibble(k = c(10, 15, 20)) %>%
  mutate(panel      = map(k, fuzzy_panel),
         n_proteins = map_int(panel, length),
         hard_n     = map_int(k, ~ length(intersect(head(rank_ttest, .x),
                                                    head(rank_rf, .x)))))

q3c_panels                       # fuzzy vs hard panel size for each k
q3c_panels$panel[[1]]            # k = 10 fuzzy panel
setdiff(q3c_panels$panel[[1]], proteins_sstar_train)   # added vs hard (Q3a)


# Q3c panels (Rule A): k = 10/15/20 -> 5/9/16 proteins (hard: 4/5/6)
#   fuzzy k = 10 adds PTN = same panel as Q3b hard k = 15
#   fuzzy panels grow ~2-3x faster than hard; k = 20 (16 proteins) risks overfitting


## STEP 2: FIT + SCORE EACH FUZZY PANEL
########################################

q3c <- q3c_panels %>%
  mutate(metrics = map(panel, fit_score)) %>%
  select(k, n_proteins, metrics) %>%
  unnest(metrics) %>%
  mutate(variant = paste0('Q3c: fuzzy k = ', k))

# comparison table
bind_rows(bench$metrics, q3a, q3b, q3c) %>%
  select(variant, n_proteins, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

# Q3c results (Rule A: top k of one list AND top 2k of the other):
#   fuzzy k = 10 (5 proteins, + PTN): acc 0.839, AUROC 0.871 = same as Q3b hard k = 15
#   fuzzy k = 15/20 (9/16 proteins): acc 0.613, AUROC 0.788/0.783, sens 0.50/0.44
#    fuzzy helps only when it adds 1 protein; larger fuzzy panels overfit







