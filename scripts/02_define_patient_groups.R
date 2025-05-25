# 02_define_patient_groups.R

# --- 1. Load necessary libraries ---
message("Loading libraries...")
suppressPackageStartupMessages({
    library(dplyr)
    library(SummarizedExperiment)
})
message("Libraries loaded.")

# --- 2. Load filtered clinical data ---
clinical_file_path <- file.path("data", "brca_clinical_filtered.rda")
message(paste("Loading filtered clinical data from:", clinical_file_path))
if (!file.exists(clinical_file_path)) {
    stop(paste("Error: Filtered clinical data file not found at", clinical_file_path, 
               ". Please run 01_data_acquisition_preparation.R first."))
}
load(clinical_file_path) # Loads clinical_data_filtered
message("Filtered clinical data loaded.")
# Rename for clarity in this script
clinical_data <- clinical_data_filtered
rm(clinical_data_filtered) # remove the original loaded name

# --- 3. Load filtered SummarizedExperiment object ---
se_file_path <- file.path("data", "brca_se_filtered.rda")
message(paste("Loading filtered SummarizedExperiment object from:", se_file_path))
if (!file.exists(se_file_path)) {
    stop(paste("Error: Filtered SummarizedExperiment file not found at", se_file_path,
               ". Please run 01_data_acquisition_preparation.R first."))
}
load(se_file_path) # Loads brca_data_filtered
message("Filtered SummarizedExperiment object loaded.")
# Rename for clarity in this script
se_object <- brca_data_filtered
rm(brca_data_filtered) # remove the original loaded name

# --- 4. Ensure clinical data is not empty ---
message("Checking if clinical data is empty...")
if (nrow(clinical_data) == 0) {
    stop("Error: The loaded clinical data is empty. Cannot proceed with group definition.")
}
message(paste("Clinical data contains", nrow(clinical_data), "patients."))

# --- 5. Calculate the median OS_time ---
message("Calculating median OS_time...")
median_os_time <- median(clinical_data$OS_time, na.rm = TRUE)
if (is.na(median_os_time)) {
    stop("Error: Median OS_time is NA. Check OS_time column in clinical data.")
}
message(paste("Median OS_time:", round(median_os_time, 2), "days."))

# --- 6. Create survival_group column ---
message("Creating 'survival_group' column in clinical data...")
clinical_data$survival_group <- ifelse(clinical_data$OS_time <= median_os_time, 
                                       "short_survivor", 
                                       "long_survivor")
clinical_data$survival_group <- factor(clinical_data$survival_group, levels = c("short_survivor", "long_survivor"))
message("'survival_group' column created.")

# --- 7. Print summary of group assignments ---
message("Summary of survival group assignments:")
print(table(clinical_data$survival_group, useNA = "ifany"))

# --- 8. Save the updated clinical data frame ---
grouped_clinical_path <- file.path("data", "brca_clinical_grouped.rda")
message(paste("Saving updated clinical data with survival groups to:", grouped_clinical_path))
save(clinical_data, file = grouped_clinical_path)
message("Updated clinical data saved.")

# --- 9. Add survival_group to colData of SummarizedExperiment object ---
message("Adding 'survival_group' to colData of the SummarizedExperiment object...")

# Ensure patient barcodes are available in clinical_data (expected 'barcode' column from script 01)
if (!"barcode" %in% colnames(clinical_data)) {
    stop("Error: 'barcode' column missing in clinical_data. Cannot reliably match to SummarizedExperiment.")
}

# Ensure colnames(se_object) (barcodes) are present in clinical_data$barcode
if (!all(colnames(se_object) %in% clinical_data$barcode)) {
    stop("Error: Not all sample barcodes from SummarizedExperiment are present in the clinical data. Data mismatch.")
}

# Create a mapping dataframe from clinical_data
survival_group_map <- clinical_data[, c("barcode", "survival_group")]

# Get current colData from se_object
se_colData <- as.data.frame(colData(se_object))

# Add rownames to se_colData if they are not already the barcodes for easier merging
# (colnames(se_object) are the barcodes)
se_colData$barcode_temp_merge_key <- colnames(se_object)

# Merge survival_group into se_colData
se_colData_updated <- dplyr::left_join(se_colData, survival_group_map, by = c("barcode_temp_merge_key" = "barcode"))

# Check if merge was successful and all samples got a group
if (any(is.na(se_colData_updated$survival_group))) {
   message("Warning: Some samples in SummarizedExperiment did not receive a survival_group after merging.")
   print(head(se_colData_updated[is.na(se_colData_updated$survival_group),]))
}

# Ensure the order is the same as original colData before assignment
# The left_join should preserve the order of the left dataframe (se_colData) if barcodes are unique
# Then, assign the new column from the merged data.
colData(se_object)$survival_group <- factor(se_colData_updated$survival_group, levels = c("short_survivor", "long_survivor"))

# Remove the temporary merge key
colData(se_object)$barcode_temp_merge_key <- NULL

message("'survival_group' added to colData of SummarizedExperiment object.")
message("Summary of survival_group in SummarizedExperiment colData:")
print(table(colData(se_object)$survival_group, useNA = "ifany"))

# Verification: Check if barcodes align and groups match
# This is an important check
if (!identical(rownames(colData(se_object)), colnames(se_object))) {
    message("Warning: Rownames of colData(se_object) do not match colnames(se_object). This might indicate an issue.")
}
# Compare a few samples
# sample_barcodes_check <- head(colnames(se_object))
# for (bc in sample_barcodes_check) {
#    se_group <- colData(se_object)[bc, "survival_group"]
#    cl_group <- clinical_data[clinical_data$barcode == bc, "survival_group"]
#    message(paste("Barcode:", bc, "- SE group:", se_group, "- Clinical group:", cl_group))
#    if (!identical(as.character(se_group), as.character(cl_group))) {
#        message(paste("Mismatch for barcode:", bc))
#    }
# }


# --- 10. Save the updated SummarizedExperiment object ---
grouped_se_path <- file.path("data", "brca_se_grouped.rda")
message(paste("Saving updated SummarizedExperiment object with survival groups to:", grouped_se_path))
save(se_object, file = grouped_se_path)
message("Updated SummarizedExperiment object saved.")

message("\n--- Patient group definition script finished. ---")
message(paste("Updated clinical data saved to:", grouped_clinical_path))
message(paste("Updated SummarizedExperiment object saved to:", grouped_se_path))

# --- End of script ---
