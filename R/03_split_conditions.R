# Split conditions

library(here)
library(yaml)
library(data.table)
library(dplyr)

setwd(here())
config <- read_yaml("config/config.yaml")

input_file <- "data/processed/all_processed.csv"
skin_data <- fread(input_file, na.strings = c("", "NA", "NaN", "null", "NULL"))
skin_data <- as.data.frame(skin_data)

diagnosis_col <- config$data$diagnosis_col
malignant_labels <- config$labels$malignant

if (!diagnosis_col %in% names(skin_data)) {
  stop("Diagnosis column '", diagnosis_col, "' not found in dataset.")
}

unique_diagnoses <- unique(skin_data[[diagnosis_col]])
unique_diagnoses <- unique_diagnoses[!is.na(unique_diagnoses)]

skin_data$condition <- ifelse(skin_data[[diagnosis_col]] %in% malignant_labels, "malignant", "benign")

malignant_data <- skin_data[skin_data$condition == "malignant", ]
benign_data <- skin_data[skin_data$condition == "benign", ]

malignant_data <- malignant_data[, setdiff(names(malignant_data), "condition")]
benign_data <- benign_data[, setdiff(names(benign_data), "condition")]

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

write.csv(malignant_data, "data/processed/malignant.csv", row.names = FALSE)
write.csv(benign_data, "data/processed/benign.csv", row.names = FALSE)
