# 04_survival_analysis_DEGs.R

# --- 1. Load necessary libraries ---
message("Loading libraries...")
suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(SummarizedExperiment)
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
})
message("Libraries loaded.")

# --- 2. Create directory for Kaplan-Meier plots ---
km_plot_dir <- file.path("results", "KM_plots")
if (!dir.exists(km_plot_dir)) {
  dir.create(km_plot_dir, recursive = TRUE)
  message(paste("Created directory:", km_plot_dir))
} else {
  message(paste("Directory already exists:", km_plot_dir))
}

# --- 3. Load significant DEGs table ---
significant_degs_path <- file.path("results", "significant_degs.csv")
message(paste("Loading significant DEGs from:", significant_degs_path))
if (!file.exists(significant_degs_path)) {
    stop(paste("Error: Significant DEGs file not found at", significant_degs_path,
               ". Please run 03_differential_expression.R first."))
}
significant_degs_df <- read.csv(significant_degs_path, row.names = 1) # Assuming first column is gene names
if (nrow(significant_degs_df) == 0) {
    stop(paste("Error: No significant DEGs found in", significant_degs_path,
               ". Cannot proceed with survival analysis for DEGs."))
}
message(paste("Loaded", nrow(significant_degs_df), "significant DEGs."))
# Ensure padj is numeric
significant_degs_df$padj <- as.numeric(significant_degs_df$padj)
significant_degs_df$log2FoldChange <- as.numeric(significant_degs_df$log2FoldChange)


# --- 4. Load the dds object ---
dds_path <- file.path("data", "dds_object.rda")
message(paste("Loading dds object from:", dds_path))
if (!file.exists(dds_path)) {
    stop(paste("Error: dds object file not found at", dds_path,
               ". Please run 03_differential_expression.R first."))
}
load(dds_path) # Loads dds
message("dds object loaded.")

# --- 5. Extract normalized expression data ---
message("Extracting normalized expression data using VST...")
# Check if dds has enough samples for VST, otherwise use rlog
if (ncol(dds) < 30 && ncol(dds) > 1) { # rlog often preferred for smaller sample sizes
    message("Using rlog transformation as sample size < 30.")
    rld <- rlog(dds, blind = FALSE)
    norm_counts <- assay(rld)
} else if (ncol(dds) <=1) {
    stop("Error: Cannot perform VST or rlog on a DESeqDataSet with 1 or fewer samples.")
} else {
    message("Using VST transformation.")
    vst_data <- vst(dds, blind = FALSE)
    norm_counts <- assay(vst_data)
}
message("Normalized expression data extracted.")
# print(head(norm_counts[,1:5])) # Optional: view first few rows/cols

# --- 6. Extract clinical data ---
message("Extracting clinical data from dds object...")
clinical_info <- as.data.frame(colData(dds))
# Ensure OS_time and OS_status are present
if (!all(c("OS_time", "OS_status") %in% colnames(clinical_info))) {
    stop("Error: 'OS_time' and/or 'OS_status' not found in colData(dds). Check previous scripts.")
}
# Ensure they are numeric
clinical_info$OS_time <- as.numeric(as.character(clinical_info$OS_time))
clinical_info$OS_status <- as.numeric(as.character(clinical_info$OS_status))

# Filter out any samples with NA OS_time or OS_status, or OS_time <= 0 (should have been done in script 01)
valid_survival_samples <- !is.na(clinical_info$OS_time) & 
                          !is.na(clinical_info$OS_status) & 
                          clinical_info$OS_time > 0
clinical_info <- clinical_info[valid_survival_samples, ]
norm_counts <- norm_counts[, rownames(clinical_info)] # Ensure norm_counts match filtered clinical_info

if (nrow(clinical_info) == 0) {
    stop("Error: No samples remaining after filtering for valid survival data in clinical_info.")
}
message(paste("Clinical data extracted for", nrow(clinical_info), "samples with valid survival information."))


# --- 7. Select top DEGs for analysis ---
message("Selecting top DEGs for survival analysis...")
# Filter for padj < 0.05 (already done if significant_degs.csv is used correctly)
degs_for_survival <- significant_degs_df[significant_degs_df$padj < 0.05 & !is.na(significant_degs_df$padj), ]

# Separate up-regulated (in long survivors) and down-regulated (in long survivors)
up_regulated_degs <- degs_for_survival %>% 
                        filter(log2FoldChange > 0) %>% 
                        arrange(padj, desc(log2FoldChange))
down_regulated_degs <- degs_for_survival %>% 
                        filter(log2FoldChange < 0) %>% 
                        arrange(padj, log2FoldChange) # most negative LFC

# Select top 5 from each group
top_n <- 5
selected_up_degs <- head(up_regulated_degs, top_n)
selected_down_degs <- head(down_regulated_degs, top_n)

selected_degs_for_km <- rbind(selected_up_degs, selected_down_degs)
# Add gene names from rownames if they are not a column
if (!"gene_name" %in% colnames(selected_degs_for_km) && identical(rownames(selected_degs_for_km), character(0))==FALSE) {
    selected_degs_for_km$gene_name <- rownames(selected_degs_for_km)
} else if (!"gene_name" %in% colnames(selected_degs_for_km) && "X" %in% colnames(selected_degs_for_km)){
    # If gene names were in first column of CSV and read as "X"
    selected_degs_for_km$gene_name <- selected_degs_for_km$X
}


if (nrow(selected_degs_for_km) == 0) {
    message("No DEGs met the criteria for selection (padj < 0.05, and top N up/down). Skipping KM plot generation.")
} else {
    message(paste("Selected", nrow(selected_degs_for_km), "DEGs for Kaplan-Meier analysis:"))
    print(selected_degs_for_km[, c("gene_name", "log2FoldChange", "padj")])
}

# --- 8. Perform survival analysis for each selected DEG ---
survival_summary_list <- list()

if (nrow(selected_degs_for_km) > 0) {
    for (i in 1:nrow(selected_degs_for_km)) {
        gene_id <- selected_degs_for_km$gene_name[i]
        message(paste("\nProcessing gene:", gene_id, "(", i, "of", nrow(selected_degs_for_km), ")"))

        # a. Get normalized expression values
        if (!gene_id %in% rownames(norm_counts)) {
            message(paste("Warning: Gene ID", gene_id, "not found in normalized counts matrix. Skipping."))
            next
        }
        gene_expression <- norm_counts[gene_id, ]

        # b. Stratify patients by median expression
        median_expr <- median(gene_expression, na.rm = TRUE)
        if (is.na(median_expr)) {
            message(paste("Warning: Median expression for gene", gene_id, "is NA. Skipping."))
            next
        }
        expression_group <- ifelse(gene_expression > median_expr, "high_expression", "low_expression")
        expression_group <- factor(expression_group, levels = c("low_expression", "high_expression"))

        # Combine with clinical data for analysis
        # Ensure barcodes/sample IDs match between expression data and clinical_info
        # rownames(clinical_info) should be the sample identifiers matching colnames(norm_counts)
        if(!identical(names(gene_expression), rownames(clinical_info))){
            message("Warning: Sample ID mismatch between expression data and clinical info for gene ", gene_id, ". Attempting to align.")
            common_samples <- intersect(names(gene_expression), rownames(clinical_info))
            if(length(common_samples) == 0) {
                message("Error: No common samples after alignment for gene ", gene_id, ". Skipping.")
                next
            }
            gene_expression_aligned <- gene_expression[common_samples]
            clinical_info_aligned <- clinical_info[common_samples, ]
            expression_group_aligned <- expression_group[match(common_samples, names(gene_expression))]
        } else {
            gene_expression_aligned <- gene_expression
            clinical_info_aligned <- clinical_info
            expression_group_aligned <- expression_group
        }
        
        if(length(unique(expression_group_aligned)) < 2){
            message(paste("Warning: Gene", gene_id, "has less than 2 expression groups after stratification (e.g., all high or all low due to ties at median or zero variance). Skipping KM plot."))
            next
        }

        analysis_data <- data.frame(
            OS_time = clinical_info_aligned$OS_time,
            OS_status = clinical_info_aligned$OS_status,
            expression_group = expression_group_aligned
        )
        
        # c. Create Surv object
        surv_obj <- Surv(time = analysis_data$OS_time, event = analysis_data$OS_status)

        # d. Perform survival analysis
        fit <- survfit(surv_obj ~ expression_group, data = analysis_data)

        # e. Generate Kaplan-Meier plot
        # Sanitize gene_id for filename
        sanitized_gene_id <- gsub("[^A-Za-z0-9_.-]", "_", gene_id)
        plot_title <- paste("Kaplan-Meier Plot for", gene_id)
        
        km_plot <- tryCatch({
             ggsurvplot(
                fit,
                data = analysis_data,
                title = plot_title,
                pval = TRUE, # Show p-value from log-rank test
                pval.method = TRUE, # Show method for p-value calculation
                legend.title = "Expression Group",
                legend.labs = c("Low Expression", "High Expression"),
                risk.table = TRUE,
                risk.table.col = "strata",
                risk.table.y.text.col = TRUE,
                risk.table.y.text = FALSE,
                xlab = "Time (days)",
                ylab = "Overall Survival Probability",
                palette = c("#E7B800", "#2E9FDF"), # Example colors
                ggtheme = theme_minimal()
            )
        }, error = function(e) {
            message(paste("Error generating KM plot for", gene_id, ":", e$message))
            return(NULL)
        })

        if (!is.null(km_plot)) {
            # f. Save the plot
            plot_filename <- file.path(km_plot_dir, paste0("KM_plot_", sanitized_gene_id, ".pdf"))
            tryCatch({
                pdf(plot_filename, width = 8, height = 7)
                print(km_plot, newpage = FALSE) # newpage=FALSE for ggsurvplot object
                dev.off()
                message(paste("Saved KM plot to:", plot_filename))
            }, error = function(e_save) {
                message(paste("Error saving KM plot for", gene_id, "to PDF:", e_save$message))
                # Attempt PNG as fallback
                plot_filename_png <- file.path(km_plot_dir, paste0("KM_plot_", sanitized_gene_id, ".png"))
                tryCatch({
                    ggsave(plot_filename_png, print(km_plot), width = 8, height = 7)
                    message(paste("Saved KM plot to (PNG fallback):", plot_filename_png))
                }, error = function(e_save_png){
                    message(paste("Error saving KM plot for", gene_id, "to PNG:", e_save_png$message))
                })
            })

            # Store p-value for summary
            surv_diff <- survdiff(surv_obj ~ expression_group, data = analysis_data)
            p_value <- pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE)
            survival_summary_list[[gene_id]] <- data.frame(
                gene = gene_id,
                log2FoldChange = selected_degs_for_km[selected_degs_for_km$gene_name == gene_id, "log2FoldChange"],
                padj_DEA = selected_degs_for_km[selected_degs_for_km$gene_name == gene_id, "padj"],
                log_rank_p_value = p_value,
                median_expr_cutoff = median_expr
            )
        }
    }
}

# --- 9. Create and save summary table ---
if (length(survival_summary_list) > 0) {
    survival_summary_df <- do.call(rbind, survival_summary_list)
    rownames(survival_summary_df) <- NULL # Clean rownames
    summary_table_path <- file.path("results", "survival_analysis_summary_DEGs.csv")
    message(paste("\nSaving survival analysis summary table to:", summary_table_path))
    write.csv(survival_summary_df, file = summary_table_path, row.names = FALSE)
    message("Survival analysis summary table saved.")
    print(survival_summary_df)
} else {
    message("\nNo survival analyses were successfully performed or no DEGs selected, so no summary table generated.")
}

message("\n--- Survival analysis script for DEGs finished. ---")
message(paste("Kaplan-Meier plots saved in:", km_plot_dir))
if (length(survival_summary_list) > 0) {
    message(paste("Survival summary table saved to:", file.path("results", "survival_analysis_summary_DEGs.csv")))
}
# --- End of script ---
