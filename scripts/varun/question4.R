# Question 4: find a simpler panel with comparable accuracy, or a better
# panel with improved accuracy, benchmarked against the in-class pipeline.
#
# Benchmark (from scripts/chris/chris-inclasspipelinenotes.R, verified
# reproducible): panel = DERM, RELT, MRC2, IgD, Cadherin-5 (5 proteins)
# sens 0.812, spec 0.733, acc 0.774, AUROC 0.883
#
# In that panel's logistic regression, only DERM (p=0.020) and IgD
# (p=0.004) were statistically significant -- RELT, MRC2, and Cadherin-5
# were not. That's the lead for a simpler panel: try DERM + IgD alone.

library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)
library(glmnet)
load('data/biomarker-clean.RData')

benchmark <- readRDS('results/benchmark.rds')
proteins_sstar <- benchmark$sstar
proteins_sstar

## Recreate the exact same train/test split as the benchmark
## (same data, same seed, same order of operations) so the comparison
## is apples-to-apples -- only the protein panel differs.
biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

## SIMPLER PANEL: just the 2 statistically significant proteins
#################################################################
simple_proteins <- c('DERM', 'IgD')

train_simple <- training(biomarker_split) %>% select(class, all_of(simple_proteins))
test_simple  <- testing(biomarker_split) %>% select(class, all_of(simple_proteins))

fit_simple <- glm(class ~ ., data = train_simple, family = 'binomial')
summary(fit_simple)

simple_results <- test_simple %>%
  add_predictions(fit_simple, type = 'response') %>%
  mutate(truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred, event_level = 'second') %>%
  mutate(variant = 'simple panel (DERM + IgD)', n_proteins = length(simple_proteins))

simple_results

## COMPARISON: benchmark (5 proteins) vs. simple panel (2 proteins)
#####################################################################
comparison <- bind_rows(
  benchmark$metrics %>% select(.metric, .estimate, variant),
  simple_results %>% select(.metric, .estimate, variant)
) %>%
  pivot_wider(names_from = variant, values_from = .estimate)

comparison

saveRDS(simple_results, 'results/question4_simple_panel.rds')

## LASSO: let the model pick the panel from ALL 1,317 proteins at once,
## instead of pre-filtering with t-tests + random forest first.
## Same 80/20 split as the benchmark -- initial_split only depends on
## nrow() and the seed, not which columns are present, so using the
## full protein set here still reproduces the identical train/test rows.
#####################################################################
biomarker_full <- biomarker_clean %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group, -ados)   # ados excluded from predictors, same as in-class pipeline

set.seed(101422)
split_full <- biomarker_full %>% initial_split(prop = 0.8)

train_full <- training(split_full)
test_full  <- testing(split_full)

x_train <- train_full %>% select(-class) %>% as.matrix()
y_train <- train_full$class
x_test  <- test_full %>% select(-class) %>% as.matrix()

# 5-fold cross-validation to choose the penalty strength (lambda)
set.seed(101422)
cv_lasso <- cv.glmnet(x_train, y_train, family = 'binomial', alpha = 1, nfolds = 5)

# lambda.min: best cross-validated fit (more proteins kept)
coefs_min <- coef(cv_lasso, s = 'lambda.min')
selected_min <- setdiff(rownames(coefs_min)[which(coefs_min[, 1] != 0)], '(Intercept)')
selected_min
length(selected_min)

pred_min <- predict(cv_lasso, newx = x_test, s = 'lambda.min', type = 'response')[, 1]

lasso_min_results <- test_full %>%
  mutate(pred = pred_min,
         truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred, event_level = 'second') %>%
  mutate(variant = 'LASSO (lambda.min)', n_proteins = length(selected_min))

lasso_min_results

# lambda.1se: more conservative, fewer proteins kept
coefs_1se <- coef(cv_lasso, s = 'lambda.1se')
selected_1se <- setdiff(rownames(coefs_1se)[which(coefs_1se[, 1] != 0)], '(Intercept)')
selected_1se
length(selected_1se)

pred_1se <- predict(cv_lasso, newx = x_test, s = 'lambda.1se', type = 'response')[, 1]

lasso_1se_results <- test_full %>%
  mutate(pred = pred_1se,
         truth = factor(class, levels = c(FALSE, TRUE)),
         est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
  class_metrics(truth = truth, estimate = est, pred, event_level = 'second') %>%
  mutate(variant = 'LASSO (lambda.1se)', n_proteins = length(selected_1se))

lasso_1se_results

## FULL COMPARISON: benchmark vs. simple panel vs. both LASSO variants
#########################################################################
full_comparison <- bind_rows(
  benchmark$metrics %>% select(.metric, .estimate, variant),
  simple_results %>% select(.metric, .estimate, variant),
  lasso_min_results %>% select(.metric, .estimate, variant),
  lasso_1se_results %>% select(.metric, .estimate, variant)
) %>%
  pivot_wider(names_from = variant, values_from = .estimate)

full_comparison

saveRDS(list(lasso_min = lasso_min_results, selected_min = selected_min,
             lasso_1se = lasso_1se_results, selected_1se = selected_1se),
        'results/question4_lasso.rds')

## ---------------------------------------------------------------------
## FURTHER EXPLORATION: everything above is based on ONE 31-person test
## split. With only 31 test subjects, flipping a single classification
## moves accuracy by ~3 percentage points -- too noisy to trust small
## differences between panels. Redo the comparison with repeated 10-fold
## cross-validation so every subject gets used as a test case multiple
## times, giving a much more stable estimate of "which panel is best."
## ---------------------------------------------------------------------

eval_panel_cv <- function(proteins, data, v = 10, repeats = 10, seed = 101422) {
  set.seed(seed)
  folds <- vfold_cv(data %>% select(class, all_of(proteins)), v = v, repeats = repeats)

  fold_metrics <- map_dfr(folds$splits, function(split) {
    tr <- training(split)
    te <- testing(split)
    fit <- glm(class ~ ., data = tr, family = 'binomial')
    te %>%
      mutate(pred = predict(fit, newdata = te, type = 'response'),
             truth = factor(class, levels = c(FALSE, TRUE)),
             est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
      class_metrics(truth = truth, estimate = est, pred, event_level = 'second')
  })

  fold_metrics %>%
    group_by(.metric) %>%
    summarise(mean = mean(.estimate), sd = sd(.estimate), .groups = 'drop')
}

union_proteins <- union(benchmark$s1, benchmark$s2)   # 10 + 10 minus overlap
length(union_proteins)

panels <- list(
  'benchmark (5: t-test ∩ RF)'      = proteins_sstar,
  'simple (2: DERM + IgD)'          = c('DERM', 'IgD'),
  'add RELT back (3)'               = c('DERM', 'IgD', 'RELT'),
  'union (15: t-test ∪ RF)'         = union_proteins
)

cv_comparison <- imap_dfr(panels, function(proteins, name) {
  eval_panel_cv(proteins, biomarker_full) %>% mutate(panel = name, n_proteins = length(proteins))
})

cv_comparison %>%
  select(panel, n_proteins, .metric, mean, sd) %>%
  arrange(.metric, desc(mean)) %>%
  print(n = Inf)

cv_wide <- cv_comparison %>%
  select(panel, n_proteins, .metric, mean) %>%
  pivot_wider(names_from = .metric, values_from = mean)

cv_wide

saveRDS(cv_comparison, 'results/question4_cv_comparison.rds')

## Nested CV for LASSO too, for a fair comparison -- refit LASSO fresh
## in each outer fold (not just once on the original 80/20 split), so
## its number is on the same honest footing as the panels above.
#####################################################################
set.seed(101422)
folds_full <- vfold_cv(biomarker_full, v = 10)

lasso_cv_metrics <- map_dfr(folds_full$splits, function(split) {
  tr <- training(split)
  te <- testing(split)
  xtr <- tr %>% select(-class) %>% as.matrix()
  ytr <- tr$class
  xte <- te %>% select(-class) %>% as.matrix()

  fit_cv <- cv.glmnet(xtr, ytr, family = 'binomial', alpha = 1, nfolds = 5)
  pred <- predict(fit_cv, newx = xte, s = 'lambda.min', type = 'response')[, 1]

  te %>%
    mutate(pred = pred,
           truth = factor(class, levels = c(FALSE, TRUE)),
           est   = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
    class_metrics(truth = truth, estimate = est, pred, event_level = 'second')
})

lasso_cv_summary <- lasso_cv_metrics %>%
  group_by(.metric) %>%
  summarise(mean = mean(.estimate), sd = sd(.estimate), .groups = 'drop') %>%
  mutate(panel = 'LASSO (nested CV, lambda.min)', n_proteins = NA)

lasso_cv_summary

## FINAL ROBUST COMPARISON -- everything on equal footing (CV-based)
#####################################################################
final_comparison <- bind_rows(cv_comparison, lasso_cv_summary) %>%
  select(panel, n_proteins, .metric, mean) %>%
  pivot_wider(names_from = .metric, values_from = mean) %>%
  arrange(desc(accuracy))

final_comparison

saveRDS(final_comparison, 'results/question4_final_comparison.rds')
