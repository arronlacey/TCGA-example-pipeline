# TCGA Differential Gene Expression and Survival Analysis Pipeline

## Description

This pipeline performs an end-to-end analysis to identify differentially expressed genes (DEGs) between patient groups with different survival outcomes using data from The Cancer Genome Atlas (TCGA). It then conducts survival analysis on the identified DEGs. The primary cancer type used for demonstration is TCGA-BRCA (Breast Invasive Carcinoma), but the pipeline can be adapted for other cancer types.

## Prerequisites

1.  **R:** You need to have R installed on your system. You can download it from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).
2.  **R Packages:** The pipeline requires several R packages. These can be installed by running the first script of the pipeline (`scripts/00_install_packages.R`). The required packages include:
    *   `BiocManager`
    *   `TCGAbiolinks`
    *   `SummarizedExperiment`
    *   `DESeq2`
    *   `dplyr`
    *   `survival`
    *   `survminer`
    *   `ggplot2`
    *   `pheatmap`
    *   `ggrepel`
    *   `R.utils`

## Directory Structure

-   `scripts/`: Contains all the R scripts for the pipeline.
-   `data/`: Stores downloaded and processed data files (e.g., TCGA expression data, clinical information, R objects). This directory is created by the scripts.
-   `results/`: Stores output files from the analysis, including differential expression tables, lists of significant DEGs, and plots (volcano plot, heatmap). This directory is created by the scripts.
-   `results/KM_plots/`: Stores Kaplan-Meier survival plots for individual DEGs. This directory is created by the scripts.
-   `GDCdata_BRCA/` (or similar, depending on project): Default directory created by `TCGAbiolinks` for storing raw downloaded TCGA files.

## Pipeline Workflow

The R scripts are designed to be run in sequence. You can execute them using an R IDE or from the command line (e.g., `Rscript scripts/script_name.R`).

1.  **`scripts/00_install_packages.R`**:
    *   Installs all necessary R packages. Run this script first.
2.  **`scripts/01_data_acquisition_preparation.R`**:
    *   Downloads RNA-Seq (HTSeq counts) and clinical data for the specified TCGA project (default: TCGA-BRCA).
    *   Preprocesses the data: matches samples, filters based on survival information, and performs initial gene filtering.
    *   Outputs: `brca_rnaseq_prepared.rda`, `brca_se_filtered.rda`, `brca_clinical_filtered.rda` in the `data/` directory.
    *   *Note: Data download can take a significant amount of time and disk space.*
3.  **`scripts/02_define_patient_groups.R`**:
    *   Loads processed clinical data.
    *   Defines patient survival groups (e.g., 'short_survivor', 'long_survivor') based on median Overall Survival (OS) time.
    *   Outputs: `brca_clinical_grouped.rda`, `brca_se_grouped.rda` (updated `SummarizedExperiment` with group info) in `data/`.
4.  **`scripts/03_differential_expression.R`**:
    *   Performs differential gene expression analysis using DESeq2 between the defined survival groups.
    *   Outputs:
        *   `results/diff_expr_results_all.csv`: Full table of DESeq2 results.
        *   `results/significant_degs.csv`: Table of significantly differentially expressed genes (padj < 0.05, |log2FC| > 1).
        *   `data/dds_object.rda`, `data/diff_expr_results_ordered.rda`: DESeq2 objects for further analysis.
5.  **`scripts/04_survival_analysis_DEGs.R`**:
    *   Conducts Kaplan-Meier survival analysis for selected top DEGs.
    *   For each selected DEG, it stratifies patients into high/low expression groups and compares survival.
    *   Outputs:
        *   Individual Kaplan-Meier plots (e.g., `KM_plot_GENENAME.pdf`) in `results/KM_plots/`.
        *   `results/survival_analysis_summary_DEGs.csv`: Summary table of log-rank p-values for analyzed DEGs.
6.  **`scripts/05_visualization_reporting.R`**:
    *   Generates summary visualizations.
    *   Outputs:
        *   `results/volcano_plot_DEGs.pdf`: Volcano plot of differential expression results.
        *   `results/heatmap_top_DEGs.pdf`: Heatmap of top significant DEGs.

## Customization

*   **Cancer Type:** To analyze a different TCGA cancer cohort, change the `project_id` variable in `scripts/01_data_acquisition_preparation.R`.
*   **Thresholds:** Significance thresholds for p-values and log2 fold changes, as well as the number of genes for plots, can be adjusted within the respective R scripts (`03_differential_expression.R`, `04_survival_analysis_DEGs.R`, `05_visualization_reporting.R`).

---
*This pipeline provides a template for bioinformatics analysis and may require adjustments based on specific research questions and data characteristics.*
