# Data Download Guide for SCINA Workflow

## Quick Reference

**Best Dataset:** AD03501 (SEA-AD Middle Temporal Gyrus Atlas)
- **Size:** ~378,000 cells
- **Format:** .qsave (Seurat object)
- **Source:** ssREAD database

---

## Method 1: Manual Download from ssREAD (Recommended)

### Step-by-Step Instructions

1. **Navigate to ssREAD Downloads**
   - Open browser: https://bmblx.bmi.osumc.edu/ssread/
   - Click "Downloads" in top menu

2. **Browse Datasets**
   - Look for "sc/snRNA-seq Data Downloads" section
   - Filter by:
     - Species: **Human**
     - Assay: **snRNA-seq**
     - Region: **Middle Temporal Gyrus (MTG)** or **Prefrontal Cortex (PFC)**

3. **Find Dataset AD03501**
   - Project ID: **AD03501**
   - Description: SEA-AD Middle Temporal Gyrus Atlas
   - Cells: ~378,211
   - File: AD03501.qsave

4. **Download**
   - Click the download link
   - Save to: `Fa25-Project6-AD-Transcriptomics/data/ssread_datasets/`
   - Expected file size: **~2-3 GB** (if smaller, download failed)

5. **Verify Download**
   ```r
   # In R:
   file.size("data/ssread_datasets/AD03501.qsave") / 1e9  # Should be ~2-3 GB
   ```

6. **Load in R**
   ```r
   library(qs)
   seurat_obj <- qread("data/ssread_datasets/AD03501.qsave")
   seurat_obj  # Check object structure
   ```

---

## Method 2: Alternative Datasets

### Option A: Smaller Test Dataset

If AD03501 is too large for your computer:

**AD001 Series** (ROSMAP cohort)
- Individual samples: AD00101, AD00102, etc.
- ~5,000-10,000 cells each
- Prefrontal cortex
- Easier to process

**Download:** Same steps as above, but select AD00101.qsave instead

### Option B: Mathys et al. from Synapse

**Direct download from original source:**

1. **Create Synapse Account**
   - Visit: https://www.synapse.org/
   - Click "Register"
   - Verify email

2. **Navigate to Dataset**
   - Go to: https://www.synapse.org/#!Synapse:syn18485175
   - Read data use agreement

3. **Download via Web Interface**
   - Click "Files" tab
   - Find processed data files (.rds or .h5ad)
   - Click download

4. **Download via R** (Programmatic)
   ```r
   # Install Synapse client
   install.packages("synapser",
                    repos = c("http://ran.synapse.org",
                              "http://cran.fhcrc.org"))

   # Login
   library(synapser)
   synLogin(email = "your_email@example.com",
            password = "your_password")

   # Download
   entity <- synGet("syn18485175",
                    downloadLocation = "data/ssread_datasets/")

   # Files will be in data/ssread_datasets/syn18485175/
   ```

5. **Convert to Seurat** (if needed)
   ```r
   # If data is in .h5ad format (AnnData)
   library(Seurat)
   library(SeuratDisk)

   Convert("data/mathys.h5ad", dest = "h5seurat")
   seurat_obj <- LoadH5Seurat("data/mathys.h5seurat")

   # Save as RDS for faster loading
   saveRDS(seurat_obj, "data/ssread_datasets/mathys_2019_seurat.rds")
   ```

---

## Method 3: GEO Download (Lau et al. 2020)

**Advantage:** No account needed, well-documented

1. **Visit GEO**
   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827

2. **Download Supplementary Files**
   - Scroll to bottom: "Download family" section
   - Click: "GSE157827_RAW.tar"
   - Unpack: `tar -xvf GSE157827_RAW.tar`

3. **Load into Seurat**
   ```r
   # Files are usually in Matrix Market format (.mtx)
   library(Seurat)

   # Read 10x-style data
   lau_data <- Read10X(data.dir = "data/GSE157827/")

   # Create Seurat object
   seurat_obj <- CreateSeuratObject(counts = lau_data,
                                     project = "Lau2020")

   # Save
   saveRDS(seurat_obj, "data/ssread_datasets/lau_2020_seurat.rds")
   ```

---

## Recommended Dataset Comparison

| Dataset | Source | Cells | Pros | Cons |
|---------|--------|-------|------|------|
| **AD03501** | ssREAD | 378K | Pre-processed, comprehensive | Very large (slow) |
| **AD001** | ssREAD | 10K | Small, fast | Limited power |
| **Mathys** | Synapse | 80K | Gold standard, well-annotated | Requires account |
| **Lau** | GEO | 170K | Large, no account | Requires processing |

**Our Recommendation:** Start with **Mathys from Synapse** (good balance of size/quality)

---

## Data Structure Guide

### What's in an ssREAD .qsave file?

```r
# Load file
seurat_obj <- qread("AD03501.qsave")

# Check structure
str(seurat_obj, max.level = 2)

# Key components:
# @assays$RNA@data       - Normalized expression (log-scale)
# @assays$RNA@scale.data - Scaled expression (for PCA)
# @meta.data             - Cell metadata (condition, sex, etc.)
# @reductions$pca        - PCA coordinates
# @reductions$umap       - UMAP coordinates
# @meta.data$celltype    - Original cell type labels

# IMPORTANT: ssREAD objects may LACK raw counts (@assays$RNA@counts)
# This is due to data sharing restrictions
# You can still use normalized data for SCINA and visualization
```

### Required Metadata Columns

For full analysis, check if these exist:

```r
colnames(seurat_obj@meta.data)

# Ideal columns:
# - "condition"    : AD vs Control
# - "sex"          : M vs F
# - "age"          : Age at death
# - "braak_stage"  : AD severity (0-6)
# - "apoe_genotype": APOE variant (ε3/ε3, ε3/ε4, ε4/ε4)
# - "celltype"     : Pre-existing annotations (for validation)
```

---

## Troubleshooting Downloads

### Problem: Download link doesn't work

**Cause:** ssREAD is a dynamic web app, direct links may not work

**Solution:** Must use web browser and manually click download buttons

### Problem: Downloaded file is tiny (<10 MB)

**Cause:** Downloaded an error page (HTML) instead of data

**Solution:**
1. Delete the small file
2. Try again with web browser (not wget/curl)
3. Ensure JavaScript is enabled
4. Try different browser (Chrome, Firefox)

### Problem: Can't load .qsave file in R

**Error:** `Error: 'qread' not found`

**Solution:**
```r
install.packages("qs")
library(qs)
seurat_obj <- qread("file.qsave")
```

**Error:** `Error: corrupt data frame`

**Solution:** File download was interrupted
```r
# Re-download the file
# Check file size matches expected (~2-3 GB for AD03501)
file.size("AD03501.qsave") / 1e9
```

### Problem: Synapse login fails

**Error:** `Invalid username or password`

**Solutions:**
1. Check email verification (Synapse sends confirmation email)
2. Reset password on Synapse website
3. Use API key instead of password:
   ```r
   # Get API key from: https://www.synapse.org/#!PersonalAccessTokens:
   synLogin(authToken = "your_api_token_here")
   ```

### Problem: Out of memory loading large file

**Error:** `cannot allocate vector of size X GB`

**Solutions:**

1. **Increase R memory** (Windows)
   ```r
   memory.limit(size = 16000)  # 16 GB
   ```

2. **Subset on load** (if using h5ad)
   ```r
   library(anndata)
   adata <- read_h5ad("file.h5ad", backed = "r")  # Read-only mode
   # Subset before converting to Seurat
   ```

3. **Use smaller dataset** (AD001 instead of AD03501)

4. **Process on server/cloud**
   - Google Colab (free, 12 GB RAM)
   - RStudio Cloud

---

## Quick Validation

After downloading any dataset, run this validation:

```r
# Load data
seurat_obj <- readRDS("your_file.rds")  # or qread for .qsave

# 1. Check size
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")

# 2. Check metadata
cat("Metadata columns:\n")
print(colnames(seurat_obj@meta.data))

# 3. Check our genes of interest
genes_of_interest <- c("ACE", "GRN", "APOE", "APP", "APH1B")
genes_present <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
cat("\nGenes present:", paste(genes_present, collapse = ", "), "\n")

# 4. Quick visualization
library(Seurat)
if ("umap" %in% names(seurat_obj@reductions)) {
  DimPlot(seurat_obj, reduction = "umap")
} else {
  cat("No UMAP found - will be created during workflow\n")
}

# 5. Check if ready for SCINA
cat("\nData slot (for SCINA):", class(seurat_obj@assays$RNA@data), "\n")
cat("Range:", range(seurat_obj@assays$RNA@data), "\n")
cat("Expected: Matrix, range ~0-10\n")
```

**Expected output:**
```
Cells: 80000
Genes: 25000
Metadata columns:
 [1] "orig.ident"    "nCount_RNA"    "nFeature_RNA"  "condition"
 [5] "sex"           "age"           "braak_stage"

Genes present: GRN, APOE, APP, APH1B

Data slot (for SCINA): dgCMatrix
Range: 0 9.2
Expected: Matrix, range ~0-10
✓ Ready for SCINA!
```

---

## Summary: Fastest Path

**For immediate start:**

1. Download **Mathys et al. from Synapse** (syn18485175)
2. Or use **demo mode** in SCINA_workflow.R (creates test data automatically)
3. Replace with real data once downloaded

**Script handles missing data gracefully:**
- If no file found, creates demo dataset
- Runs full workflow on demo (for testing)
- Replace demo with real data later and re-run

**No blockers!** You can start running the analysis immediately.

---

## Need Help?

**Slack:** `#project-6-rna-seq`

**Common issues:**
- Can't download: Try different browser
- Can't load: Check file size, reinstall `qs` package
- Out of memory: Use smaller dataset first

---

**Ready to proceed?** Go back to `SCINA_workflow.R` and run!
