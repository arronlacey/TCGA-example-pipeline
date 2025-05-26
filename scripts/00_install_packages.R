# Check for and install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Define Bioconductor packages
bioc_packages <- c("TCGAbiolinks", "DESeq2", "SummarizedExperiment")

# Define CRAN packages
cran_packages <- c("survival", "survminer", "ggplot2", "pheatmap", "R.utils", "rmarkdown", "dplyr", "ggrepel", "knitr")

# Install Bioconductor packages
message("Checking and installing Bioconductor packages...")
for (pkg in bioc_packages) {
    message(paste("Checking for", pkg, "..."))
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        if (requireNamespace(pkg, quietly = TRUE)) {
            message(paste(pkg, "installed successfully."))
        } else {
            message(paste("Failed to install", pkg))
        }
    } else {
        message(paste(pkg, "is already installed."))
    }
}

# Install CRAN packages
message("\nChecking and installing CRAN packages...")
for (pkg in cran_packages) {
    message(paste("Checking for", pkg, "..."))
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        install.packages(pkg)
        if (requireNamespace(pkg, quietly = TRUE)) {
            message(paste(pkg, "installed successfully."))
        } else {
            message(paste("Failed to install", pkg))
        }
    } else {
        message(paste(pkg, "is already installed."))
    }
}

message("\nPackage installation check complete.")
