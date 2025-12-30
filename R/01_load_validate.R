# Load and validate

library(here)
library(yaml)
library(data.table)
library(janitor)
library(dplyr)

setwd(here())
config <- read_yaml("config/config.yaml")
dir.create("results", showWarnings = FALSE, recursive = TRUE)

validation_log <- character()

input_file <- config$data$input_csv
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

skin_data <- fread(input_file, na.strings = c("", "NA", "NaN", "null", "NULL"))
skin_data <- clean_names(skin_data)
skin_data <- as.data.frame(skin_data)

validation_log <- c(validation_log, paste("Dataset loaded from:", input_file))
validation_log <- c(validation_log, paste("Load timestamp:", Sys.time()))

n_obs <- nrow(skin_data)
n_vars <- ncol(skin_data)

var_types <- sapply(skin_data, function(x) class(x)[1])
var_type_counts <- table(var_types)

missing_counts <- colSums(is.na(skin_data))
missing_pct <- round(100 * missing_counts / n_obs, 2)
missing_summary <- data.frame(
  variable = names(missing_counts),
  n_missing = missing_counts,
  pct_missing = missing_pct
) %>%
  arrange(desc(n_missing))

validation_log <- c(validation_log, paste("n =", n_obs, "rows"))
validation_log <- c(validation_log, paste("p =", n_vars, "variables"))
validation_log <- c(validation_log, "\nVariable types:")
for (type in names(var_type_counts)) {
  validation_log <- c(validation_log, paste("  ", type, ":", var_type_counts[type]))
}
validation_log <- c(validation_log, "\nMissing values:")
validation_log <- c(validation_log, paste("  Total variables:", n_vars))
validation_log <- c(validation_log, paste("  Variables with missing:", sum(missing_counts > 0)))
validation_log <- c(validation_log, paste("  Variables complete:", sum(missing_counts == 0)))

diagnosis_col <- config$data$diagnosis_col
if (!diagnosis_col %in% names(skin_data)) {
  stop("Diagnosis column '", diagnosis_col, "' not found in dataset.")
}

diagnosis_values <- unique(skin_data[[diagnosis_col]])
diagnosis_values <- diagnosis_values[!is.na(diagnosis_values)]

malignant_labels <- config$labels$malignant
malignant_found <- diagnosis_values[diagnosis_values %in% malignant_labels]
benign_found <- diagnosis_values[!diagnosis_values %in% malignant_labels]

unknown_labels <- setdiff(diagnosis_values, c(malignant_labels, "ACK", "NEV", "SCC", "BCC", "MEL"))
if (length(unknown_labels) > 0) {
  validation_log <- c(validation_log, paste("\nWARNING: Unknown diagnosis labels:", 
                                            paste(unknown_labels, collapse = ", ")))
}

validation_log <- c(validation_log, paste("\nDiagnosis column:", diagnosis_col))
validation_log <- c(validation_log, paste("  Total unique values:", length(diagnosis_values)))
validation_log <- c(validation_log, paste("  Malignant labels:", paste(malignant_labels, collapse = ", ")))
validation_log <- c(validation_log, paste("  Malignant found in data:", paste(malignant_found, collapse = ", ")))
validation_log <- c(validation_log, paste("  Benign/other found:", paste(benign_found, collapse = ", ")))

if (length(config$variables$include) > 0) {
  missing_vars <- setdiff(config$variables$include, names(skin_data))
  if (length(missing_vars) > 0) {
    stop("Variables specified in config not found in dataset: ", 
         paste(missing_vars, collapse = ", "))
  }
  validation_log <- c(validation_log, paste("\nSpecified variables:", length(config$variables$include), "found"))
} else {
  id_cols <- config$data$id_col
  exclude_cols <- c(id_cols, diagnosis_col)
  candidate_vars <- setdiff(names(skin_data), exclude_cols)
  validation_log <- c(validation_log, paste("\nInferred candidate variables:", length(candidate_vars)))
  validation_log <- c(validation_log, paste("  Excluded (IDs + diagnosis):", length(exclude_cols)))
}

id_cols <- config$data$id_col
exclude_cols <- c(id_cols, diagnosis_col)
analysis_vars <- setdiff(names(skin_data), exclude_cols)

type_inference <- data.frame(
  variable = analysis_vars,
  current_type = var_types[analysis_vars],
  n_unique = sapply(skin_data[analysis_vars], function(x) length(unique(x[!is.na(x)]))),
  n_missing = missing_counts[analysis_vars],
  stringsAsFactors = FALSE
) %>%
  mutate(
    inferred_type = case_when(
      current_type %in% c("integer", "numeric") & n_unique > 10 ~ "continuous",
      current_type %in% c("integer", "numeric") & n_unique <= 10 ~ "ordinal",
      current_type == "logical" ~ "binary",
      current_type == "character" & n_unique <= 20 ~ "categorical",
      current_type == "character" & n_unique > 20 ~ "categorical",
      TRUE ~ "unknown"
    )
  )

validation_log <- c(validation_log, "\nType inference (preview):")
for (i in 1:nrow(type_inference)) {
  validation_log <- c(validation_log, 
    paste("  ", type_inference$variable[i], ":", 
          type_inference$current_type[i], "->", type_inference$inferred_type[i]))
}

writeLines(validation_log, "results/validation_summary.txt")
