library(tidyverse)

biomarker_data <- read.csv('data/biomarker-raw.csv')

head(biomarker_data)

dim(biomarker_data)

biomarker_data[1:5, 1:8]
unique(biomarker_data[[3]])[1:20] 

# Removing the first row of protein abbreviations
biomarker_data <- biomarker_data[-1, ]

# Converting protein columns to numeric types
biomarker_data[, 3:ncol(biomarker_data)] <-
  lapply(biomarker_data[, 3:ncol(biomarker_data)], as.numeric)

# Checking classes
sapply(biomarker_data[, 1:8], class)

# Should be 155 1320
dim(biomarker_data)

# Histograms of 4 proteins
par(mfrow = c(2, 2))

hist(biomarker_data[[3]],
     main = names(biomarker_data)[3],
     xlab = "Raw protein levels")

hist(biomarker_data[[4]],
     main = names(biomarker_data)[4],
     xlab = "Raw protein levels")

hist(biomarker_data[[5]],
     main = names(biomarker_data)[5],
     xlab = "Raw protein levels")

hist(biomarker_data[[6]],
     main = names(biomarker_data)[6],
     xlab = "Raw protein levels")



