# 05_visualization_reporting.R

# --- 1. Load necessary libraries ---
message("Loading libraries...")
suppressPackageStartupMessages({
    library(ggplot2)
    library(pheatmap)
    library(dplyr)
    library(SummarizedExperiment)
    library(DESeq2)
    library(ggrepel) # For non-overlapping labels on volcano plot
})
message("Libraries loaded.")

# --- 2. Create results directory (if it doesn't exist) ---
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  message(paste("Created directory:", results_dir))
} else {
  message(paste("Directory already exists:", results_dir))
}

# --- 3. Load full DESeq2 results ---
results_path <- file.path("data", "diff_expr_results_ordered.rda")
message(paste("Loading full DESeq2 results from:", results_path))
if (!file.exists(results_path)) {
    stop(paste("Error: DESeq2 results file not found at", results_path,
               ". Please run 03_differential_expression.R first."))
}
load(results_path) # Loads res_ordered
all_results_df <- as.data.frame(res_ordered)
all_results_df$gene <- rownames(res_ordered) # Add gene names as a column
message("Full DESeq2 results loaded.")
# print(head(all_results_df)) # Optional: view first few rows

# --- 4. Load the dds object ---
dds_path <- file.path("data", "dds_object.rda")
message(paste("Loading dds object from:", dds_path))
if (!file.exists(dds_path)) {
    stop(paste("Error: dds object file not found at", dds_path,
               ". Please run 03_differential_expression.R first."))
}
load(dds_path) # Loads dds
message("dds object loaded.")

# --- 5. Generate and Save Volcano Plot ---
message("Generating Volcano Plot...")
# a. Create a column to indicate significance
padj_threshold <- 0.05
log2fc_threshold <- 1.0
all_results_df <- all_results_df %>%
    mutate(
        is_significant = ifelse(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold, "Significant", "Not Significant"),
        is_significant = factor(is_significant, levels = c("Significant", "Not Significant")) # Control color order
    ) %>%
    filter(!is.na(padj)) # Remove rows with NA padj for plotting

# Select top genes to label (e.g., top 10 significant by padj, then by LFC magnitude)
genes_to_label_df <- all_results_df %>%
    filter(is_significant == "Significant") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(10)

volcano_plot <- ggplot(all_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = is_significant), alpha = 0.6, size = 1.5) +
    geom_text_repel(
        data = genes_to_label_df,
        aes(label = gene),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.5, "lines"),
        segment.color = 'grey50',
        max.overlaps = Inf # Allow more overlaps if necessary, or adjust
    ) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(
        title = "Volcano Plot of Differentially Expressed Genes",
        subtitle = paste("Comparison: long_survivor vs. short_survivor"),
        x = "log2(Fold Change)",
        y = "-log10(Adjusted P-value)",
        color = "Significance"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    ) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "blue")

volcano_plot_path <- file.path(results_dir, "volcano_plot_DEGs.pdf")
tryCatch({
    ggsave(volcano_plot_path, volcano_plot, width = 8, height = 7, device = "pdf")
    message(paste("Volcano plot saved to:", volcano_plot_path))
}, error = function(e) {
    message(paste("Error saving volcano plot:", e$message))
})


# --- 6. Generate and Save Heatmap of Top DEGs ---
message("Generating Heatmap of Top DEGs...")
# a. Get normalized expression data
# Using VST as it's generally good for visualization and robust to sample size differences
message("Applying VST transformation for heatmap...")
if (ncol(dds) < 30 && ncol(dds) > 1) {
    message("Using rlog transformation for heatmap as sample size < 30.")
    rld <- rlog(dds, blind = FALSE) # blind=FALSE as we have a design
    norm_counts_heatmap <- assay(rld)
} else if (ncol(dds) <=1) {
    stop("Error: Cannot perform VST or rlog on a DESeqDataSet with 1 or fewer samples.")
} else {
    message("Using VST transformation for heatmap.")
    vst_data <- vst(dds, blind = FALSE)
    norm_counts_heatmap <- assay(vst_data)
}


# b. Identify top DEGs
significant_degs_df <- all_results_df %>%
    filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold & !is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange)))

# Select top N DEGs (e.g., up to 50, or fewer if not that many significant)
num_top_degs <- min(50, nrow(significant_degs_df))
if (num_top_degs == 0) {
    message("No significant DEGs found to plot heatmap. Skipping heatmap generation.")
} else {
    message(paste("Selecting top", num_top_degs, "DEGs for heatmap."))
    top_degs_for_heatmap <- head(significant_degs_df, n = num_top_degs)
    
    # Ensure gene names are valid and exist in norm_counts_heatmap
    valid_genes <- top_degs_for_heatmap$gene[top_degs_for_heatmap$gene %in% rownames(norm_counts_heatmap)]
    if(length(valid_genes) == 0) {
        message("None of the selected top DEGs are found in the normalized counts matrix. Skipping heatmap.")
    } else {
        if(length(valid_genes) < num_top_degs) {
            message(paste("Warning: Only", length(valid_genes), "out of", num_top_degs, "selected DEGs were found in the expression matrix."))
        }
        top_degs_matrix <- norm_counts_heatmap[valid_genes, , drop = FALSE]

        # c. Prepare column annotations
        if (!"survival_group" %in% colnames(colData(dds))) {
            message("Warning: 'survival_group' not found in colData(dds). Heatmap annotation will be missing this.")
            annotation_col_df <- data.frame(row.names = colnames(dds)) # Empty annotation
        } else {
            annotation_col_df <- as.data.frame(colData(dds)[, "survival_group", drop = FALSE])
            rownames(annotation_col_df) <- colnames(dds) # Ensure rownames match sample names
        }

        # d. Use pheatmap
        heatmap_plot_path <- file.path(results_dir, "heatmap_top_DEGs.pdf")
        
        # Adjust font size based on number of genes
        fontsize_row_val <- if(num_top_degs > 30) 6 else if(num_top_degs > 20) 8 else 10

        tryCatch({
            pheatmap(
                top_degs_matrix,
                annotation_col = annotation_col_df,
                scale = "row", # Scale genes to have mean 0 and SD 1
                cluster_rows = TRUE,
                cluster_cols = TRUE, # Cluster samples based on expression of these genes
                show_colnames = FALSE, # Usually too many samples to show names
                show_rownames = TRUE,
                fontsize_row = fontsize_row_val,
                main = paste("Heatmap of Top", length(valid_genes), "Differentially Expressed Genes"),
                filename = heatmap_plot_path,
                width = 8, height = 10
            )
            message(paste("Heatmap saved to:", heatmap_plot_path))
        }, error = function(e) {
            message(paste("Error generating or saving heatmap:", e$message))
        })
    }
}

message("\n--- Visualization and reporting script finished. ---")
message(paste("Plots saved in:", results_dir))

# --- End of script ---
