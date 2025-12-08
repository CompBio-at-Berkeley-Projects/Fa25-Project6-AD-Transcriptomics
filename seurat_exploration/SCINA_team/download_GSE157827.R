###############################################################################
# Robust end-to-end pipeline for GSE157827 (final patched script)
# - Validates 10x folders
# - Extracts or downloads the full tarball if needed
# - Re-downloads individual corrupted sample files via curl
# - Loads only valid samples into Seurat, merges, preprocesses, runs PCA/UMAP/clustering
# - Saves merged Seurat object and QC plots
#
# Usage: source("download_GSE157827_final.R")
###############################################################################

## ------------ Configuration ------------
setwd("~/Documents/GitHub/Fa25-Project6-AD-Transcriptomics")  # change if needed

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Paths
raw_parent <- "data/GSE157827"
extract_dir <- file.path(raw_parent, "GSE157827_RAW")
tar_file <- file.path(raw_parent, "GSE157827_RAW.tar")
out_rds <- file.path(raw_parent, "GSE157827_processed", "GSE157827_merged_seurat.rds")
fig_dir <- file.path(raw_parent, "figures_GSE157827_QC")
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Expected sample IDs (based on this dataset)
expected_samples <- c("AD1","AD2","AD4","AD5","AD6","AD8","AD9","AD10","AD13","AD19","AD20")

message("Starting pipeline...")

## ------------ Helper functions ------------
validate_sample <- function(sample_dir) {
  bf <- file.path(sample_dir, "barcodes.tsv.gz")
  ff <- file.path(sample_dir, "features.tsv.gz")
  mf <- file.path(sample_dir, "matrix.mtx.gz")
  exists_bar <- file.exists(bf)
  exists_feat <- file.exists(ff)
  exists_mtx <- file.exists(mf)
  mtx_read <- FALSE
  dims <- ""
  emsg <- ""
  if (exists_mtx) {
    tryCatch({
      mm <- Matrix::readMM(mf)
      mtx_read <- TRUE
      dims <- paste0(nrow(mm), " x ", ncol(mm))
    }, error = function(e) {
      mtx_read <<- FALSE
      emsg <<- e$message
    })
  }
  data.frame(
    sample = basename(sample_dir),
    barcodes_ok = exists_bar,
    features_ok = exists_feat,
    matrix_ok = exists_mtx,
    mtx_readable = mtx_read,
    mtx_dim = dims,
    error_message = emsg,
    stringsAsFactors = FALSE
  )
}

download_sample_files_via_curl <- function(gsm_id, sample_name, dest_dir) {
  # gsm_id like "GSM4775598"
  # sample_name like "AD20"
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  base_path <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/",
                      substr(gsm_id, 1, 7), "nnn/", gsm_id, "/suppl/")
  files <- c(
    barcodes = paste0(gsm_id, "_", sample_name, "_barcodes.tsv.gz"),
    features = paste0(gsm_id, "_", sample_name, "_features.tsv.gz"),
    matrix   = paste0(gsm_id, "_", sample_name, "_matrix.mtx.gz")
  )
  for (k in names(files)) {
    url <- paste0(base_path, files[[k]])
    dest <- file.path(dest_dir, ifelse(k=="matrix","matrix.mtx.gz", paste0(k, ".tsv.gz")))
    # use curl with user-agent and follow redirects
    cmd <- sprintf("curl -L -A 'Mozilla/5.0' '%s' -o '%s' --retry 3 --retry-delay 5", url, dest)
    message("Downloading ", sample_name, " ", k, " from: ", url)
    rc <- system(cmd)
    if (rc != 0 || !file.exists(dest)) {
      message("Warning: download failed for ", sample_name, " ", k, " (rc=", rc, ")")
    } else {
      message("Downloaded: ", dest)
    }
  }
}

ensure_tar_extracted <- function() {
  # If extract_dir has expected sample folders, done.
  if (dir.exists(extract_dir) && length(list.dirs(extract_dir, recursive = FALSE))>0) {
    return(TRUE)
  }
  # If tar file exists, extract it using system tar (more reliable)
  if (file.exists(tar_file)) {
    message("Extracting existing tar file: ", tar_file)
    system(sprintf("tar -xvf '%s' -C '%s'", tar_file, raw_parent))
    return(TRUE)
  }
  # Otherwise, attempt to download tar via HTTPS using curl (robust)
  message("Tar file not found. Attempting to download full tar (1.2 GB) via curl...")
  url_tar <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157827/suppl/GSE157827_RAW.tar"
  dir.create(raw_parent, recursive = TRUE, showWarnings = FALSE)
  cmd <- sprintf("curl -L -A 'Mozilla/5.0' '%s' -o '%s' --retry 3 --retry-delay 10", url_tar, tar_file)
  rc <- system(cmd)
  if (rc != 0) {
    message("Warning: tar download failed (rc=", rc, "). Please download manually from GEO and place at: ", tar_file)
    return(FALSE)
  }
  message("Download complete; extracting tar...")
  system(sprintf("tar -xvf '%s' -C '%s'", tar_file, raw_parent))
  return(TRUE)
}

## ------------ Ensure extracted files exist (extract or download tar if necessary) ------------
if (!dir.exists(extract_dir) || length(list.dirs(extract_dir, recursive = FALSE)) == 0) {
  ok <- ensure_tar_extracted()
  if (!ok) {
    message("Proceeding with whatever folders are present under ", extract_dir)
  }
}

message("Listing sample directories under: ", extract_dir)
present_dirs <- list.dirs(extract_dir, recursive = FALSE, full.names = FALSE)
message("Found: ", paste(present_dirs, collapse = ", "))

# Use only the intersection of expected_samples and present_dirs; but allow any folder that looks valid
candidate_dirs <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
# Filter to folders that contain barcodes/features/matrix
valid_folder_mask <- sapply(candidate_dirs, function(d) {
  all(c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz") %in% list.files(d))
})
candidate_dirs <- candidate_dirs[valid_folder_mask]

message("Candidate 10x folders: ", paste(basename(candidate_dirs), collapse = ", "))

if (length(candidate_dirs) == 0) {
  message("No complete 10x folders found yet under ", extract_dir, ". Attempting to re-extract or re-download tar.")
  ok <- ensure_tar_extracted()
  candidate_dirs <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
  valid_folder_mask <- sapply(candidate_dirs, function(d) {
    all(c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz") %in% list.files(d))
  })
  candidate_dirs <- candidate_dirs[valid_folder_mask]
  if (length(candidate_dirs) == 0) {
    stop("No valid sample folders found after extraction. Please ensure raw 10x folders are present under: ", extract_dir)
  }
}

## ------------ Validate samples (fast check and readMM check) ------------
qc_list <- lapply(candidate_dirs, validate_sample)
qc_df <- bind_rows(qc_list)
print(qc_df)

good_samples <- qc_df %>% filter(mtx_readable==TRUE) %>% pull(sample)
bad_samples <- qc_df %>% filter(mtx_readable==FALSE) %>% pull(sample)

message("Good samples: ", paste(good_samples, collapse = ", "))
message("Bad/Corrupt samples: ", paste(bad_samples, collapse = ", "))

## ------------ If any bad samples found, attempt targeted re-download (only for AD20 known GSM) ------------
# Map sample->GSM id if known (only AD20 here; expand if you have others)
sample_to_gsm <- list(AD20 = "GSM4775598")  # add more as needed

for (s in bad_samples) {
  if (s %in% names(sample_to_gsm)) {
    gsm <- sample_to_gsm[[s]]
    message("Attempting to re-download files for ", s, " (", gsm, ")")
    download_sample_files_via_curl(gsm, s, file.path(extract_dir, s))
    # Re-validate
    new_qc <- validate_sample(file.path(extract_dir, s))
    print(new_qc)
    if (isTRUE(new_qc$mtx_readable)) {
      message("Re-download fixed sample: ", s)
      # update qc_df and lists
      qc_df <- qc_df %>% mutate(
        mtx_readable = ifelse(sample==s, TRUE, mtx_readable),
        mtx_dim = ifelse(sample==s, new_qc$mtx_dim, mtx_dim),
        error_message = ifelse(sample==s, new_qc$error_message, error_message)
      )
      good_samples <- unique(c(good_samples, s))
      bad_samples <- setdiff(bad_samples, s)
    } else {
      message("Re-download DID NOT fix ", s, " — you will need to manually download files from GEO.")
    }
  } else {
    message("No automatic GSM mapping for ", s, ". Please download sample files manually from GEO.")
  }
}

# Final lists after possible re-downloads
message("Final good samples: ", paste(good_samples, collapse = ", "))
message("Final bad samples (will be skipped): ", paste(bad_samples, collapse = ", "))

if (length(good_samples) == 0) stop("No readable samples available. Fix data files and re-run.")

## ------------ Load good samples into Seurat ------------
seurat_list <- list()
for (s in good_samples) {
  sample_dir <- file.path(extract_dir, s)
  message("Loading sample: ", s)
  # Guarded Read10X
  tryCatch({
    counts <- Read10X(sample_dir)
    obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
    obj$sample <- s
    obj$condition <- ifelse(grepl("^AD", s), "AD", "Control")
    seurat_list[[s]] <- obj
    message("Loaded ", s, " (cells=", ncol(obj), ", genes=", nrow(obj), ")")
  }, error = function(e) {
    message("ERROR loading sample ", s, ": ", e$message)
  })
}

if (length(seurat_list) == 0) stop("No Seurat objects created. Abort.")

## ------------ Merge safely ------------
message("Merging ", length(seurat_list), " Seurat objects...")
if (length(seurat_list) == 1) {
  merged <- seurat_list[[1]]
} else {
  merged <- Reduce(function(a,b) merge(a, y=b), seurat_list)
}
message("Merged object: cells=", ncol(merged), " genes=", nrow(merged))

## ------------ QC metrics and filtering ------------
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

qc_stats <- data.frame(
  Metric = c("Total cells", "Total genes", "Median genes/cell", "Median UMIs/cell", "Median mito %"),
  Value = c(ncol(merged), nrow(merged), median(merged$nFeature_RNA), median(merged$nCount_RNA), median(merged$percent.mt))
)
print(qc_stats)

# Filter thresholds (customize if you want)
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
message("After QC filter: cells=", ncol(merged), " genes=", nrow(merged))

## ------------ Preprocessing ------------
message("Normalizing...")
merged <- NormalizeData(merged, verbose = FALSE)

message("Finding variable features...")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

message("Scaling data...")
merged <- ScaleData(merged, verbose = FALSE)

message("Running PCA (50 PCs)...")
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

message("Running UMAP (dims 1:30)...")
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

message("Finding neighbors and clusters...")
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.4, verbose = FALSE)

## ------------ Plots & Save ------------
p1 <- DimPlot(merged, group.by = "sample") + ggtitle("GSE157827: by sample")
p2 <- DimPlot(merged, label = TRUE) + ggtitle("GSE157827: clusters")

ggsave(file.path(fig_dir, "umap_by_sample.pdf"), p1, width = 8, height = 6)
ggsave(file.path(fig_dir, "umap_clusters.pdf"), p2, width = 8, height = 6)

saveRDS(merged, out_rds)
message("Saved merged Seurat object: ", out_rds)
message("Saved QC plots to: ", fig_dir)

############################################################
## 12. Save as qsave (qs format) — NEW ADDITION
############################################################

qs_dir <- file.path(raw_parent, "qs_processed")
dir.create(qs_dir, recursive = TRUE, showWarnings = FALSE)

qs_file <- file.path(qs_dir, "GSE157827_merged_seurat.qs")

message("Saving qsave (qs) file...")
qsave(merged, qs_file, preset = "high")
message("✔ qsave complete: ", qs_file)

############################################################

## ------------ Final summary ------------
message("Pipeline complete.")
message("Good samples: ", paste(good_samples, collapse = ", "))
message("Skipped bad samples: ", paste(bad_samples, collapse = ", "))
message("Merged cells: ", ncol(merged), ", genes: ", nrow(merged))

###############################################################################
# End of script
###############################################################################
