################################################################################
#
# AUTOMATED DOWNLOAD SCRIPT FOR GSE157827 (Lau et al. 2020)
#
# This script downloads and prepares the complete GSE157827 dataset
# (169,496 nuclei from human AD prefrontal cortex) for SCINA analysis.
#
# ADVANTAGES over ssREAD:
# - Complete count matrices (not redacted)
# - Direct download (no manual clicking)
# - Standard 10x format (easy to load)
# - All target genes present (APOE, APP, GRN, ACE, APH1B)
#
# USAGE:
#   source("seurat_exploration/SCINA_team/download_GSE157827.R")
#
# TIME: ~15-30 minutes (depending on internet speed)
#
################################################################################

cat("============================================================\n")
cat("GSE157827 Download and Preparation Script\n")
cat("Lau et al. 2020 - AD Prefrontal Cortex snRNA-seq\n")
cat("============================================================\n\n")

# Set working directory
setwd("~/Documents/GitHub/Fa25-Project6-AD-Transcriptomics")

# ==============================================================================
# SECTION 1: SETUP
# ==============================================================================

cat("SECTION 1: Setup and Dependencies\n")
cat("---------------------------------\n\n")

# Create directories
cat("Creating directories...\n")
dir.create("data/GSE157827", showWarnings = FALSE, recursive = TRUE)
dir.create("data/GSE157827_processed", showWarnings = FALSE, recursive = TRUE)

# Install/load required packages
required_packages <- c("Seurat", "dplyr", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg)
  }
}

library(Seurat)
library(dplyr)
library(ggplot2)

cat("Setup complete!\n\n")

# ==============================================================================
# SECTION 2: DOWNLOAD DATA
# ==============================================================================

cat("SECTION 2: Downloading GSE157827 from GEO\n")
cat("------------------------------------------\n\n")

# File paths
tar_file <- "data/GSE157827/GSE157827_RAW.tar"
extract_dir <- "data/GSE157827/GSE157827_RAW"

# Check if already downloaded
if (file.exists(tar_file)) {
  cat("TAR file already exists. Skipping download.\n")
  cat(sprintf("File size: %.2f GB\n", file.size(tar_file) / 1e9))
} else {
  cat("Downloading GSE157827_RAW.tar (1.2 GB)...\n")
  cat("This may take 10-20 minutes depending on your connection.\n\n")

  # FTP download URL
  ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157827/suppl/GSE157827_RAW.tar"

  # Try download with error handling
  tryCatch({
    download.file(
      url = ftp_url,
      destfile = tar_file,
      method = "auto",  # auto-detect best method
      mode = "wb",       # binary mode
      quiet = FALSE
    )
    cat("\nDownload complete!\n")
  }, error = function(e) {
    cat("\n\nERROR: Automatic download failed.\n")
    cat("Error message:", e$message, "\n\n")
    cat("MANUAL DOWNLOAD INSTRUCTIONS:\n")
    cat("1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827\n")
    cat("2. Click 'GSE157827_RAW.tar (1.2 Gb)' under Supplementary file\n")
    cat("3. Save to:", tar_file, "\n")
    cat("4. Re-run this script\n\n")
    stop("Download failed - follow manual instructions above")
  })
}

# ==============================================================================
# SECTION 3: EXTRACT TAR FILE
# ==============================================================================

cat("\n")
cat("SECTION 3: Extracting TAR file\n")
cat("-------------------------------\n\n")

if (dir.exists(extract_dir) && length(list.files(extract_dir)) > 0) {
  cat("Files already extracted. Skipping extraction.\n")
} else {
  cat("Extracting TAR file...\n")
  cat("This will create ~63 files (21 samples x 3 files each)\n\n")

  # Extract TAR
  untar(tarfile = tar_file, exdir = extract_dir)

  # Check extraction
  extracted_files <- list.files(extract_dir)
  cat(sprintf("Extracted %d files\n", length(extracted_files)))

  if (length(extracted_files) == 0) {
    stop("ERROR: Extraction failed. Check if TAR file is complete.")
  }
}

cat("\nExtraction complete!\n\n")

# ==============================================================================
# SECTION 4: ORGANIZE FILES BY SAMPLE
# ==============================================================================

cat("SECTION 4: Organizing files by sample\n")
cat("--------------------------------------\n\n")

# List all files
all_files <- list.files(extract_dir, full.names = TRUE)

# Parse sample names
# File format: GSM4775579_AD1_barcodes.tsv.gz
sample_ids <- unique(sapply(strsplit(basename(all_files), "_"), function(x) x[2]))
sample_ids <- sample_ids[sample_ids != ""]  # Remove empty strings

cat(sprintf("Found %d samples:\n", length(sample_ids)))

# Separate AD and Control
ad_samples <- grep("^AD", sample_ids, value = TRUE)
nc_samples <- grep("^NC", sample_ids, value = TRUE)

cat(sprintf("  - AD samples: %d (%s)\n", length(ad_samples), paste(head(ad_samples, 3), collapse = ", ")))
cat(sprintf("  - Control samples: %d (%s)\n", length(nc_samples), paste(head(nc_samples, 3), collapse = ", ")))

# Create sample directories
cat("\nCreating sample directories...\n")
for (sample in sample_ids) {
  sample_dir <- file.path(extract_dir, sample)
  dir.create(sample_dir, showWarnings = FALSE)

  # Move files for this sample
  sample_files <- grep(paste0("_", sample, "_"), all_files, value = TRUE)

  for (file in sample_files) {
    file_type <- tail(strsplit(basename(file), "_")[[1]], 1)  # barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz
    new_name <- file.path(sample_dir, file_type)

    if (!file.exists(new_name)) {
      file.rename(file, new_name)
    }
  }
}

cat("Files organized!\n\n")

# ==============================================================================
# SECTION 5: LOAD AND MERGE ALL SAMPLES
# ==============================================================================

cat("SECTION 5: Loading data into Seurat\n")
cat("------------------------------------\n\n")

# Check if merged object already exists
merged_rds <- "data/GSE157827_processed/GSE157827_merged_seurat.rds"

if (file.exists(merged_rds)) {
  cat("Merged Seurat object already exists. Loading...\n")
  seurat_obj <- readRDS(merged_rds)
  cat(sprintf("Loaded: %d cells, %d genes\n\n", ncol(seurat_obj), nrow(seurat_obj)))

} else {
  cat("Loading samples into Seurat...\n")
  cat("This may take 5-10 minutes for all 21 samples.\n\n")

  seurat_list <- list()

  for (sample in sample_ids) {
    sample_dir <- file.path(extract_dir, sample)

    cat(sprintf("Loading %s...", sample))

    # Check if all required files exist
    required_files <- c("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")
    files_exist <- sapply(required_files, function(f) {
      file.exists(file.path(sample_dir, f))
    })

    if (all(files_exist)) {
      tryCatch({
        # Read 10x data
        data <- Read10X(data.dir = sample_dir)

        # Create Seurat object
        seurat_sample <- CreateSeuratObject(
          counts = data,
          project = "Lau2020",
          min.cells = 3,
          min.features = 200
        )

        # Add metadata
        seurat_sample$sample <- sample
        seurat_sample$condition <- ifelse(grepl("^AD", sample), "AD", "Control")

        seurat_list[[sample]] <- seurat_sample

        cat(sprintf(" %d cells\n", ncol(seurat_sample)))

      }, error = function(e) {
        cat(sprintf(" ERROR: %s\n", e$message))
      })
    } else {
      cat(" SKIPPED (missing files)\n")
    }
  }

  # Merge all samples
  cat("\nMerging all samples...\n")
  seurat_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "GSE157827"
  )

  cat(sprintf("Merged dataset: %d cells, %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))

  # Add aggregate metadata
  seurat_obj$dataset <- "GSE157827_Lau2020"
  seurat_obj$brain_region <- "Prefrontal_Cortex"
  seurat_obj$technology <- "10x_Genomics"

  # Save merged object
  cat("\nSaving merged Seurat object...\n")
  saveRDS(seurat_obj, merged_rds)
  cat(sprintf("Saved: %s\n\n", merged_rds))
}

# ==============================================================================
# SECTION 6: VERIFY GENE PRESENCE
# ==============================================================================

cat("SECTION 6: Verifying target genes\n")
cat("----------------------------------\n\n")

genes_of_interest <- c("APOE", "APP", "GRN", "ACE", "APH1B")

cat("Checking for target genes:\n")
for (gene in genes_of_interest) {
  present <- gene %in% rownames(seurat_obj)
  status <- if (present) "✓ FOUND" else "✗ MISSING"
  cat(sprintf("  %-10s %s\n", gene, status))
}

genes_present <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
genes_missing <- genes_of_interest[!genes_of_interest %in% rownames(seurat_obj)]

cat(sprintf("\nResult: %d/%d genes present\n", length(genes_present), length(genes_of_interest)))

if (length(genes_missing) > 0) {
  cat("\nSearching for missing genes with alternative names...\n")
  for (gene in genes_missing) {
    matches <- grep(gene, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) {
      cat(sprintf("  Possible matches for %s: %s\n", gene, paste(matches, collapse = ", ")))
    }
  }
}

cat("\n")

# ==============================================================================
# SECTION 7: QUICK QC SUMMARY
# ==============================================================================

cat("SECTION 7: Quality Control Summary\n")
cat("-----------------------------------\n\n")

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC metrics
qc_stats <- data.frame(
  Metric = c("Total cells", "Total genes", "Median genes/cell",
             "Median UMIs/cell", "Median mito %",
             "AD samples", "Control samples"),
  Value = c(
    ncol(seurat_obj),
    nrow(seurat_obj),
    round(median(seurat_obj$nFeature_RNA)),
    round(median(seurat_obj$nCount_RNA)),
    round(median(seurat_obj$percent.mt), 2),
    sum(seurat_obj$condition == "AD"),
    sum(seurat_obj$condition == "Control")
  )
)

print(qc_stats)

# Sample breakdown
cat("\n\nSample breakdown:\n")
sample_counts <- table(seurat_obj$sample, seurat_obj$condition)
print(sample_counts)

cat("\n")

# ==============================================================================
# SECTION 8: BASIC PREPROCESSING
# ==============================================================================

cat("SECTION 8: Basic Preprocessing\n")
cat("------------------------------\n\n")

# Check if preprocessing already done
if (!"umap" %in% names(seurat_obj@reductions)) {
  cat("Running basic preprocessing...\n")
  cat("(Full preprocessing will be done in SCINA_workflow.R)\n\n")

  # Normalize
  cat("  - Normalizing...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Find variable features
  cat("  - Finding variable features...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                      nfeatures = 2000, verbose = FALSE)

  # Scale
  cat("  - Scaling data...\n")
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  # PCA
  cat("  - Running PCA...\n")
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

  # UMAP
  cat("  - Running UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

  # Clustering
  cat("  - Clustering...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

  # Save
  saveRDS(seurat_obj, merged_rds)
  cat("\nPreprocessing complete!\n\n")

} else {
  cat("Preprocessing already done. Skipping.\n\n")
}

# ==============================================================================
# SECTION 9: GENERATE QC PLOTS
# ==============================================================================

cat("SECTION 9: Generating QC plots\n")
cat("-------------------------------\n\n")

# Create figures directory
dir.create("figures/GSE157827_QC", showWarnings = FALSE, recursive = TRUE)

# Plot 1: UMAP by condition
p1 <- DimPlot(seurat_obj, group.by = "condition") +
  ggtitle("GSE157827: AD vs Control")
ggsave("figures/GSE157827_QC/01_umap_by_condition.pdf", p1, width = 8, height = 6)
cat("Saved: figures/GSE157827_QC/01_umap_by_condition.pdf\n")

# Plot 2: UMAP by sample
p2 <- DimPlot(seurat_obj, group.by = "sample", label = FALSE) +
  ggtitle("GSE157827: All 21 Samples") +
  theme(legend.text = element_text(size = 6))
ggsave("figures/GSE157827_QC/02_umap_by_sample.pdf", p2, width = 10, height = 8)
cat("Saved: figures/GSE157827_QC/02_umap_by_sample.pdf\n")

# Plot 3: QC violin plots
p3 <- VlnPlot(seurat_obj,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition",
              ncol = 3, pt.size = 0)
ggsave("figures/GSE157827_QC/03_qc_metrics.pdf", p3, width = 12, height = 4)
cat("Saved: figures/GSE157827_QC/03_qc_metrics.pdf\n")

# Plot 4: Target genes
if (length(genes_present) > 0) {
  p4 <- FeaturePlot(seurat_obj, features = genes_present,
                    ncol = 3, cols = c("lightgrey", "red"))
  ggsave("figures/GSE157827_QC/04_target_genes.pdf", p4,
         width = 15, height = ceiling(length(genes_present)/3) * 5)
  cat("Saved: figures/GSE157827_QC/04_target_genes.pdf\n")
}

cat("\n")

# ==============================================================================
# SECTION 10: FINAL SUMMARY
# ==============================================================================

cat("============================================================\n")
cat("DOWNLOAD AND PREPARATION COMPLETE!\n")
cat("============================================================\n\n")

cat("Dataset Summary:\n")
cat(sprintf("  - Cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  - Genes: %d\n", nrow(seurat_obj)))
cat(sprintf("  - AD samples: %d\n", length(ad_samples)))
cat(sprintf("  - Control samples: %d\n", length(nc_samples)))
cat(sprintf("  - Target genes present: %d/%d\n", length(genes_present), length(genes_of_interest)))

cat("\nFiles created:\n")
cat(sprintf("  - Merged Seurat object: %s\n", merged_rds))
cat("  - QC plots: figures/GSE157827_QC/\n")

cat("\nNext Steps:\n")
cat("1. Open SCINA_workflow.R\n")
cat("2. Update line 134 to:\n")
cat('   seurat_obj <- readRDS("data/GSE157827_processed/GSE157827_merged_seurat.rds")\n')
cat("3. Run the full SCINA pipeline\n")
cat("4. Compare results with Azimuth team\n\n")

cat("Ready for SCINA analysis!\n\n")

################################################################################
# END OF SCRIPT
################################################################################
