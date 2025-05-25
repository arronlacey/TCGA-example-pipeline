# 03_differential_expression.R

# --- 1. Load necessary libraries ---
message("Loading libraries...")
suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
    library(dplyr) # For data manipulation, if needed
})
message("Libraries loaded.")

# --- 2. Load the grouped SummarizedExperiment object ---
se_grouped_path <- file.path("data", "brca_se_grouped.rda")
message(paste("Loading grouped SummarizedExperiment object from:", se_grouped_path))
if (!file.exists(se_grouped_path)) {
    stop(paste("Error: Grouped SummarizedExperiment file not found at", se_grouped_path,
               ". Please run 02_define_patient_groups.R first."))
}
load(se_grouped_path) # Loads se_object
message("Grouped SummarizedExperiment object loaded.")
message("Dimensions of loaded se_object:")
print(dim(se_object))
# print(head(colData(se_object))) # Optional: view colData

# --- 3. Ensure survival_group is a factor and set reference level ---
message("Ensuring 'survival_group' is a factor and setting reference level...")
if (!"survival_group" %in% colnames(colData(se_object))) {
    stop("Error: 'survival_group' column not found in colData of se_object. Please run 02_define_patient_groups.R.")
}
colData(se_object)$survival_group <- factor(colData(se_object)$survival_group, levels = c("short_survivor", "long_survivor"))
message("Reference level for 'survival_group' set to 'short_survivor'.")
message("Levels of survival_group: ", paste(levels(colData(se_object)$survival_group), collapse=", "))
message("Table of survival_group factor:")
print(table(colData(se_object)$survival_group))

# --- 4. Create a DESeqDataSet object ---
message("Creating DESeqDataSet object...")
# The design formula specifies how to model the counts based on the survival_group
dds <- DESeqDataSet(se_object, design = ~ survival_group)
message("DESeqDataSet object created.")
message("Dimensions of dds object:")
print(dim(dds))

# --- 5. Perform pre-filtering ---
message("Performing pre-filtering of genes with low counts...")
# Keep rows (genes) that have at least 10 reads total across all samples.
# This is a common filtering step to reduce noise and improve performance.
# This step might have been partially covered in script 01, but DESeq2's recommendation is good to follow.
keep <- rowSums(counts(dds)) >= 10
original_gene_count <- nrow(dds)
dds <- dds[keep,]
filtered_gene_count <- nrow(dds)
message(paste("Removed", original_gene_count - filtered_gene_count, "genes with total counts < 10."))
message(paste("Number of genes remaining:", filtered_gene_count))

if (nrow(dds) == 0) {
    stop("Error: No genes remaining after pre-filtering (counts < 10). Check previous filtering steps or data quality.")
}
if (ncol(dds) == 0) {
    stop("Error: No samples in dds object. This should not happen if previous steps were successful.")
}


# --- 6. Run the DESeq analysis ---
message("Running DESeq analysis (estimateSizeFactors, estimateDispersions, nbinomWaldTest)...")
# This function performs the main DESeq2 analysis steps.
# It can take some time for large datasets.
dds <- DESeq(dds)
message("DESeq analysis complete.")

# --- 7. Extract results ---
message("Extracting results for 'long_survivor' vs 'short_survivor'...")
# contrast specifies: column_name, level_for_numerator, level_for_denominator
# This will give log2FoldChange = log2(long_survivor / short_survivor)
res <- results(dds, contrast=c("survival_group", "long_survivor", "short_survivor"))
message("Results extracted.")
message("Summary of results (res):")
summary(res) # Prints summary to console

# --- 8. Order results by adjusted p-value ---
message("Ordering results by adjusted p-value (padj)...")
res_ordered <- res[order(res$padj), ]
message("Results ordered.")
# print(head(res_ordered)) # Optional: view top results

# --- 9. Save the complete ordered results table ---
message("Saving complete ordered results table...")
# Create 'results' directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
  message("Created 'results' directory.")
}
results_all_path <- file.path("results", "diff_expr_results_all.csv")
write.csv(as.data.frame(res_ordered), file=results_all_path)
message(paste("Complete results table saved to:", results_all_path))

# --- 10. Identify significant DEGs ---
padj_threshold <- 0.05
log2fc_threshold <- 1.0
message(paste("Identifying significant DEGs (padj <", padj_threshold, "and abs(log2FoldChange) >", log2fc_threshold, ")..."))
significant_degs <- subset(res_ordered, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
message(paste(nrow(significant_degs), "significant DEGs found."))
# print(head(significant_degs)) # Optional: view significant DEGs

# --- 11. Save the table of significant DEGs ---
if (nrow(significant_degs) > 0) {
    significant_degs_path <- file.path("results", "significant_degs.csv")
    message(paste("Saving table of significant DEGs to:", significant_degs_path))
    write.csv(as.data.frame(significant_degs), file=significant_degs_path)
    message("Significant DEGs table saved.")
} else {
    message("No significant DEGs found with the current thresholds. Skipping save of significant DEGs table.")
}

# --- 12. Save DESeqDataSet and ordered results objects ---
# These can be useful for later analyses, plotting, etc.
dds_path <- file.path("data", "dds_object.rda")
message(paste("Saving DESeqDataSet object (dds) to:", dds_path))
save(dds, file = dds_path)
message("dds object saved.")

res_ordered_path <- file.path("data", "diff_expr_results_ordered.rda")
message(paste("Saving ordered results object (res_ordered) to:", res_ordered_path))
save(res_ordered, file = res_ordered_path)
message("res_ordered object saved.")

message("\n--- Differential expression analysis script finished. ---")
message(paste("Full results saved to:", results_all_path))
if (nrow(significant_degs) > 0) {
    message(paste("Significant DEGs saved to:", file.path("results", "significant_degs.csv")))
}
message(paste("DESeq2 objects (dds, res_ordered) saved in '", getwd(), "/data' directory.", sep=""))

# --- End of script ---
