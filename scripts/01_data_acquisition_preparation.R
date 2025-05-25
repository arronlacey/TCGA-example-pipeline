# 01_data_acquisition_preparation.R

# --- 0. Setup ---
# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
  message("Created 'data' directory.")
} else {
  message("'data' directory already exists.")
}

# --- 1. Load necessary libraries ---
message("Loading libraries...")
suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(dplyr) # Using dplyr for easier data manipulation
})
message("Libraries loaded.")

# --- 2. Define project and cancer type ---
project_id <- "TCGA-BRCA"
message(paste("Project ID set to:", project_id))

# --- 3. Use GDCquery to find RNA-Seq HTSeq - Counts data ---
message("Querying GDC for RNA-Seq data...")
query_rnaseq <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts"
)
message("GDCquery complete. Number of results: ", nrow(getResults(query_rnaseq)))

# --- 4. Download the data using GDCdownload ---
message("Downloading data using GDCdownload...")
# Wrap GDCdownload in a tryCatch to handle potential download errors
# Small subset for testing - comment out GDCdownload for now if data is large or to avoid re-download
# query_rnaseq_subset <- query_rnaseq
# query_rnaseq_subset$results[[1]] <- query_rnaseq$results[[1]][1:20,] # download only 20 files for test
# query_rnaseq_subset$results[[1]]$cases <- query_rnaseq_subset$results[[1]]$cases[1:20]


# It is recommended to run GDCdownload and GDCprepare interactively first
# as they can take a long time and might require user input or specific proxy settings.
# For this script, we will assume it's being run in an environment where download is possible.

# Check if a prepared file already exists to avoid re-downloading if not necessary
prepared_file_path <- file.path("data", "brca_rnaseq_prepared.rda")
if (file.exists(prepared_file_path)) {
    message(paste("Prepared data file already exists at:", prepared_file_path, ". Skipping download and preparation."))
    load(prepared_file_path, verbose = TRUE) # Loads brca_data
} else {
    message("Prepared data file not found. Proceeding with download and preparation.")
    tryCatch({
        GDCdownload(query_rnaseq, method = "api", files.per.chunk = 10, directory = "GDCdata_BRCA") # Save to specific dir
    }, error = function(e_api) {
        message(paste("API GDCdownload failed:", e_api$message))
        message("Trying GDCdownload with 'client' method.")
        tryCatch({
            GDCdownload(query_rnaseq, method = "client", files.per.chunk = 10, directory = "GDCdata_BRCA")
        }, error = function(e_client) {
            message(paste("Client GDCdownload also failed:", e_client$message))
            stop("Data download failed with both API and client methods.")
        })
    })
    message("GDCdownload complete.")

    # --- 5. Prepare the data into a SummarizedExperiment object ---
    message("Preparing data using GDCprepare...")
    brca_data <- GDCprepare(query_rnaseq, 
                            save = TRUE, 
                            save.filename = prepared_file_path,
                            directory = "GDCdata_BRCA") # Read from specific dir
    message(paste("GDCprepare complete. Data saved to:", prepared_file_path))
}

message("Dimensions of prepared SummarizedExperiment object (brca_data):")
print(dim(brca_data))

# --- 6. Extract clinical data ---
message("Extracting clinical data...")
clinical_data_full <- as.data.frame(colData(brca_data))
message("Dimensions of full clinical data:")
print(dim(clinical_data_full))
# print(head(clinical_data_full)) # Optional: view first few rows

# --- 7. Filter clinical data ---
message("Filtering clinical data...")

# Select relevant columns
relevant_cols <- c("barcode", "patient", "days_to_death", "days_to_last_follow_up", "vital_status")
# Ensure all relevant columns exist
actual_relevant_cols <- relevant_cols[relevant_cols %in% colnames(clinical_data_full)]
if (length(actual_relevant_cols) < length(relevant_cols)) {
    message("Warning: Not all expected clinical columns found. Available: ", paste(colnames(clinical_data_full), collapse=", "))
    message("Using available columns: ", paste(actual_relevant_cols, collapse=", "))
}
clinical_data_subset <- clinical_data_full[, actual_relevant_cols, drop = FALSE]

# Standardize survival information
clinical_data_subset$OS_time <- NA
clinical_data_subset$OS_status <- NA

# Calculate OS_time and OS_status
# Ensure days_to_death and days_to_last_follow_up are numeric
if("days_to_death" %in% colnames(clinical_data_subset)) {
  clinical_data_subset$days_to_death <- as.numeric(as.character(clinical_data_subset$days_to_death))
}
if("days_to_last_follow_up" %in% colnames(clinical_data_subset)) {
  clinical_data_subset$days_to_last_follow_up <- as.numeric(as.character(clinical_data_subset$days_to_last_follow_up))
}


for (i in 1:nrow(clinical_data_subset)) {
    if (!is.na(clinical_data_subset$vital_status[i])) {
        if (clinical_data_subset$vital_status[i] == "Dead" || clinical_data_subset$vital_status[i] == "deceased") { # accommodate variations
            clinical_data_subset$OS_status[i] <- 1 # Deceased
            clinical_data_subset$OS_time[i] <- clinical_data_subset$days_to_death[i]
        } else if (clinical_data_subset$vital_status[i] == "Alive" || clinical_data_subset$vital_status[i] == "alive") {
            clinical_data_subset$OS_status[i] <- 0 # Alive
            clinical_data_subset$OS_time[i] <- clinical_data_subset$days_to_last_follow_up[i]
        }
    }
}
message("Standardized OS_time and OS_status.")
message("Summary of OS_time:")
print(summary(clinical_data_subset$OS_time))
message("Table of OS_status:")
print(table(clinical_data_subset$OS_status, useNA = "ifany"))


# Filter out patients with NA values for OS_time or OS_status, or OS_time <= 0
message(paste("Number of patients before NA/invalid time filtering:", nrow(clinical_data_subset)))
clinical_data_filtered <- clinical_data_subset[
    !is.na(clinical_data_subset$OS_time) &
    !is.na(clinical_data_subset$OS_status) &
    clinical_data_subset$OS_time > 0,
]
message(paste("Number of patients after NA/invalid time filtering:", nrow(clinical_data_filtered)))
if (nrow(clinical_data_filtered) == 0) {
    stop("No patients remaining after clinical data filtering. Check vital_status, days_to_death, and days_to_last_follow_up columns.")
}

# --- 8. Filter SummarizedExperiment object ---
message("Filtering SummarizedExperiment object to match filtered clinical data...")
# Ensure barcodes from clinical data are present in the SummarizedExperiment object
# The barcodes in colData(brca_data) should match those in clinical_data_filtered$barcode
common_barcodes <- intersect(colnames(brca_data), clinical_data_filtered$barcode)
message(paste("Number of common barcodes between SE and filtered clinical data:", length(common_barcodes)))

brca_data_filtered <- brca_data[, common_barcodes]
# Also filter clinical data to ensure perfect match and order
clinical_data_filtered <- clinical_data_filtered[clinical_data_filtered$barcode %in% common_barcodes, ]

message("Dimensions of brca_data after filtering for patients in clinical_data_filtered:")
print(dim(brca_data_filtered))
message("Dimensions of clinical_data_filtered after matching with brca_data_filtered:")
print(dim(clinical_data_filtered))

# --- 9. Filter genes ---
message("Filtering genes...")
counts_matrix <- assay(brca_data_filtered) # Get the counts matrix

# Remove genes with rowSums (total counts across all samples) equal to 0
genes_all_zero <- rowSums(counts_matrix) == 0
message(paste("Number of genes with rowSums = 0:", sum(genes_all_zero)))
brca_data_filtered <- brca_data_filtered[!genes_all_zero, ]
message("Dimensions of brca_data_filtered after removing all-zero genes:")
print(dim(brca_data_filtered))

# Optional: Remove genes with very low average counts (e.g., mean counts < 1)
# This is a common filtering step to reduce noise and computational burden
mean_counts <- rowMeans(assay(brca_data_filtered))
genes_low_mean_counts <- mean_counts < 1
message(paste("Number of genes with mean counts < 1:", sum(genes_low_mean_counts)))
brca_data_filtered <- brca_data_filtered[!genes_low_mean_counts, ]
message("Dimensions of brca_data_filtered after removing low-mean-count genes:")
print(dim(brca_data_filtered))

if (nrow(brca_data_filtered) == 0) {
    stop("No genes remaining after filtering. Check count data and filtering thresholds.")
}
if (ncol(brca_data_filtered) == 0) {
    stop("No samples remaining after filtering. Check clinical data and sample matching.")
}

# --- 10. Save the filtered SummarizedExperiment object ---
filtered_se_path <- file.path("data", "brca_se_filtered.rda")
message(paste("Saving filtered SummarizedExperiment object to:", filtered_se_path))
save(brca_data_filtered, file = filtered_se_path)
message("Filtered SummarizedExperiment object saved.")

# --- 11. Save the processed and matched clinical data frame ---
filtered_clinical_path <- file.path("data", "brca_clinical_filtered.rda")
message(paste("Saving processed and matched clinical data to:", filtered_clinical_path))
save(clinical_data_filtered, file = filtered_clinical_path)
message("Processed clinical data saved.")

message("\n--- Data acquisition and preparation script finished. ---")
message(paste("Final dimensions of SummarizedExperiment (brca_data_filtered):", paste(dim(brca_data_filtered), collapse = " x ")))
message(paste("Final dimensions of clinical data (clinical_data_filtered):", paste(dim(clinical_data_filtered), collapse = " x ")))
message(paste("Output files saved in '", getwd(), "/data' directory.", sep=""))

# --- End of script ---
