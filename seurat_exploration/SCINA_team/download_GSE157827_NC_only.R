###############################################################################
# Memory-efficient NC-only pipeline for GSE157827 on macOS
###############################################################################

## ------------ Configuration ------------
setwd("~/Documents/GitHub/Fa25-Project6-AD-Transcriptomics")  # adjust as needed

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(qs)

raw_parent <- "data/GSE157827"
extract_dir <- file.path(raw_parent, "GSE157827_RAW_NC_Only")

out_rds <- file.path(raw_parent, "GSE157827_processed", "GSE157827_NC_only_seurat.rds")
fig_dir <- file.path(raw_parent, "figures_GSE157827_NC_QC")
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

qs_dir <- file.path(raw_parent, "qs_processed")
dir.create(qs_dir, recursive = TRUE, showWarnings = FALSE)

message("Starting NC-only pipeline...")

## ------------ Find all matrix files ------------
matrix_files <- list.files(
  path = extract_dir,
  pattern = "_matrix\\.mtx\\.gz$",
  full.names = TRUE
)

if(length(matrix_files) == 0) stop("No matrix files found in folder: ", extract_dir)
message("Found NC matrix files:\n", paste(basename(matrix_files), collapse="\n"))

## ------------ Load and pre-filter each sample individually ------------
seurat_list <- list()

for(matrix_file in matrix_files){
  sample_id <- sub("GSM[0-9]+_(NC[0-9]+)_matrix\\.mtx\\.gz","\\1",basename(matrix_file))
  barcodes_file <- sub("_matrix\\.mtx\\.gz$","_barcodes.tsv.gz",matrix_file)
  features_file <- sub("_matrix\\.mtx\\.gz$","_features.tsv.gz",matrix_file)
  
  message("Loading sample: ", sample_id)
  
  # Create temp folder for older Seurat versions
  tmp_dir <- file.path(tempdir(), sample_id)
  dir.create(tmp_dir, showWarnings=FALSE)
  
  # Copy files to temp folder with standard names
  file.copy(matrix_file, file.path(tmp_dir,"matrix.mtx.gz"), overwrite=TRUE)
  file.copy(barcodes_file, file.path(tmp_dir,"barcodes.tsv.gz"), overwrite=TRUE)
  file.copy(features_file, file.path(tmp_dir,"features.tsv.gz"), overwrite=TRUE)
  
  # Read counts
  counts <- Read10X(data.dir = tmp_dir)
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = counts, project = sample_id,
                            min.cells = 3, min.features = 200)
  obj$sample <- sample_id
  obj$condition <- "Control"
  
  # Pre-filter cells to reduce memory
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
  
  seurat_list[[sample_id]] <- obj
  
  # Clean up
  rm(obj)
  gc()
}

good_samples <- names(seurat_list)
if(length(good_samples) == 0) stop("No valid NC Seurat objects created.")
message("Loaded NC samples: ", paste(good_samples, collapse=", "))

## ------------ Merge in memory-efficient chunks ------------
chunk1 <- seurat_list[1:ceiling(length(seurat_list)/2)]
chunk2 <- seurat_list[(ceiling(length(seurat_list)/2)+1):length(seurat_list)]

merged1 <- Reduce(function(a,b) merge(a, y=b), chunk1)
merged2 <- Reduce(function(a,b) merge(a, y=b), chunk2)

# Remove intermediate objects to free memory
rm(chunk1, chunk2, seurat_list)
gc()

# Final merge
merged <- merge(merged1, y = merged2)
rm(merged1, merged2)
gc()

message("Merged NC object: cells=", ncol(merged), " genes=", nrow(merged))

## ------------ QC & Filtering ------------
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern="^MT-")
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
message("After QC filter: cells=", ncol(merged))

## ------------ Preprocessing ------------
merged <- NormalizeData(merged, verbose=FALSE)
merged <- FindVariableFeatures(merged, selection.method="vst", nfeatures=3000)
merged <- ScaleData(merged, verbose=FALSE)
merged <- RunPCA(merged, npcs=50, verbose=FALSE)
merged <- RunUMAP(merged, dims=1:30)
merged <- FindNeighbors(merged, dims=1:30)
merged <- FindClusters(merged, resolution=0.4)

## ------------ Plots ------------
p1 <- DimPlot(merged, group.by="sample") + ggtitle("GSE157827 NC-only: by sample")
p2 <- DimPlot(merged, label=TRUE) + ggtitle("GSE157827 NC-only: clusters")

ggsave(file.path(fig_dir, "umap_by_sample_NC.pdf"), p1, width=8, height=6)
ggsave(file.path(fig_dir, "umap_clusters_NC.pdf"), p2, width=8, height=6)

## ------------ Save objects ------------
saveRDS(merged, out_rds)
message("Saved NC-only Seurat object: ", out_rds)

qs_file <- file.path(qs_dir, "GSE157827_NC_only.qs")
qsave(merged, qs_file, preset="high")
message("✔ NC-only qsave complete: ", qs_file)

###############################################################################
# End of memory-efficient NC-only pipeline
###############################################################################
