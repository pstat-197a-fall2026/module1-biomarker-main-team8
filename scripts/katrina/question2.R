library(tidyverse)
file.exists("data/biomarker-raw.csv")
var_names <- read_csv('data/biomarker-raw.csv',
                      col_names = FALSE,
                      n_max = 2,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1,
         abbreviation = V2) %>%
  na.omit()

biomarker_untrimmed <- read_csv('data/biomarker-raw.csv',
                                skip = 2,
                                col_select = -2L,
                                col_names = c('group',
                                              'empty',
                                              pull(var_names, abbreviation),
                                              'ados'),
                                na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados),
                ~ scale(log10(.x))[, 1])) %>%
  mutate(subject_id = row_number()) %>%
  select(subject_id, group, ados, everything())

View(biomarker_untrimmed)
outlier_counts <- biomarker_untrimmed %>%
  pivot_longer(
    cols = -c(subject_id, group, ados),
    names_to = "protein",
    values_to = "level"
  ) %>%
  group_by(subject_id, group) %>%
  summarise(
    n_outliers = sum(abs(level) > 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_outliers))


View(outlier_counts)
#Subject 154 has most outlier measurements 

#comparing outlier counts bw ASD and TD groups 
group_summary <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    n_subjects = n(),
    mean_outliers = mean(n_outliers),
    median_outliers = median(n_outliers),
    max_outliers = max(n_outliers),
    .groups = "drop"
  )

View(group_summary)

#Box plot 
ggplot(outlier_counts, aes(x = group, y = n_outliers, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  labs(
    title = "Number of Outlying Protein Measurements by Group",
    x = "Group",
    y = "Number of Outliers per Subject"
  ) +
  theme_minimal() +
  guides(fill = "none")

#Top subjects 
top_subjects <- outlier_counts %>%
  slice_max(n_outliers, n = 10, with_ties = FALSE)

ggplot(top_subjects,
       aes(x = reorder(factor(subject_id), n_outliers),
           y = n_outliers,
           fill = group)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 10 Subjects with the Most Outliers",
    x = "Subject ID",
    y = "Number of Outlying Protein Measurements",
    fill = "Group"
  ) +
  theme_minimal()


#After removing the outlier trimming from preprocessing, we examined the number of protein measurements more than 3 SDs from the mean for each subject. Subject 154 (TD) had the most with 157. Of the 10 subjects with the highest outlier counts, 7 belonged to TD group.
#TD also had a higher mean outlier count than ASD, but both groups had the same median of 8.5
