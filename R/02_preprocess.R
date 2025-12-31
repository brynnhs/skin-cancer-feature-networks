# Preprocessing

library(here)
library(yaml)
library(data.table)
library(dplyr)

setwd(here())
config <- read_yaml("config/config.yaml")

original_file <- "data/processed/all_processed_raw_for_phase2.csv"
if (file.exists(original_file)) {
  skin_data <- fread(original_file, na.strings = c("", "NA", "NaN", "null", "NULL"))
} else {
  skin_data <- fread(config$data$input_csv, na.strings = c("", "NA", "NaN", "null", "NULL"))
}
skin_data <- as.data.frame(skin_data)

prefilter_file <- "data/processed/all_processed_pre_complete_case.csv"
write.csv(skin_data, prefilter_file, row.names = FALSE)

id_cols <- config$data$id_col
diagnosis_col <- config$data$diagnosis_col
exclude_cols <- c(id_cols, diagnosis_col)

if (length(config$variables$include) > 0) {
  analysis_vars <- config$variables$include
} else {
  analysis_vars <- setdiff(names(skin_data), exclude_cols)
}

missing_vars <- setdiff(analysis_vars, names(skin_data))
if (length(missing_vars) > 0) {
  stop("Analysis variables not found in dataset: ", paste(missing_vars, collapse = ", "))
}

if (length(config$variables$types$continuous) > 0 || 
    length(config$variables$types$categorical) > 0 ||
    length(config$variables$types$ordinal) > 0 ||
    length(config$variables$types$binary) > 0) {
  continuous_vars <- config$variables$types$continuous
  categorical_vars <- config$variables$types$categorical
  ordinal_vars <- config$variables$types$ordinal
  binary_vars <- config$variables$types$binary
} else {
  continuous_vars <- c()
  categorical_vars <- c()
  ordinal_vars <- c()
  binary_vars <- c()
  
  for (var in analysis_vars) {
    var_class <- class(skin_data[[var]])[1]
    n_unique <- length(unique(skin_data[[var]][!is.na(skin_data[[var]])]))
    
    if (var_class == "logical") {
      binary_vars <- c(binary_vars, var)
    } else if (var_class %in% c("integer", "numeric")) {
      if (n_unique > 10) {
        continuous_vars <- c(continuous_vars, var)
      } else {
        ordinal_vars <- c(ordinal_vars, var)
      }
    } else if (var_class == "character") {
      categorical_vars <- c(categorical_vars, var)
    }
  }
}

keep_cols <- c(exclude_cols, analysis_vars)
processed_data <- skin_data[, keep_cols, drop = FALSE]

for (var in continuous_vars) {
  if (var %in% names(processed_data)) {
    processed_data[[var]] <- as.numeric(processed_data[[var]])
  }
}

for (var in categorical_vars) {
  if (var %in% names(processed_data)) {
    processed_data[[var]] <- as.factor(processed_data[[var]])
  }
}

for (var in ordinal_vars) {
  if (var %in% names(processed_data)) {
    unique_vals <- unique(processed_data[[var]][!is.na(processed_data[[var]])])
    if (is.numeric(unique_vals)) {
      ordered_levels <- sort(unique_vals)
    } else {
      ordered_levels <- sort(as.character(unique_vals))
    }
    processed_data[[var]] <- factor(processed_data[[var]], levels = ordered_levels, ordered = TRUE)
  }
}

for (var in binary_vars) {
  if (var %in% names(processed_data)) {
    processed_data[[var]] <- as.logical(processed_data[[var]])
  }
}

if (config$missing$mode == "complete_case") {
  complete_rows <- complete.cases(processed_data[, analysis_vars])
  processed_data <- processed_data[complete_rows, ]
} else if (config$missing$mode == "simple_impute") {
  for (var in analysis_vars) {
    n_missing <- sum(is.na(processed_data[[var]]))
    if (n_missing > 0) {
      if (var %in% continuous_vars) {
        median_val <- median(processed_data[[var]], na.rm = TRUE)
        processed_data[[var]][is.na(processed_data[[var]])] <- median_val
      } else if (var %in% categorical_vars || var %in% ordinal_vars) {
        mode_val <- names(sort(table(processed_data[[var]], useNA = "no"), decreasing = TRUE))[1]
        processed_data[[var]][is.na(processed_data[[var]])] <- mode_val
      } else if (var %in% binary_vars) {
        mode_val <- names(sort(table(processed_data[[var]], useNA = "no"), decreasing = TRUE))[1]
        processed_data[[var]][is.na(processed_data[[var]])] <- as.logical(mode_val)
      }
    }
  }
} else {
  stop("Unknown missing data mode: ", config$missing$mode)
}

dir.create(dirname("data/processed/all_processed.csv"), showWarnings = FALSE, recursive = TRUE)
write.csv(processed_data, "data/processed/all_processed.csv", row.names = FALSE)
