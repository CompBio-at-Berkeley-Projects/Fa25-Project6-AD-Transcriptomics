# HighRes-AD-Transcriptomics in Human Alzheimer's Disease

**Fall 2025 | Computational Biology @ Berkeley | Project 6**

---

## What Are We Actually Doing?

We're using publicly available human brain sequencing data to understand what's happening at the cellular level in Alzheimer's Disease (AD). Specifically, we're trying to figure out:

1. **Which genes are turned on or off in different brain cell types when someone has AD?**
2. **How do specific "vulnerable" cell populations (like Microglia and Astrocytes) change during disease progression?**


## The Big Picture Goal

By the end of this project, we'll have:

- **Final Scientific Report** summarizing our findings
- **Final Presentation** (think conference-style poster/talk)
- **Fully-Documented GitHub Repository** with all our code and analysis

**Official Goal Statement:**
*"To utilize human single-cell and spatial RNA sequencing data from the ssREAD/GEO/Synapse portals to characterize the cell-type-specific transcriptional profiles in Alzheimer's Disease, with a focus on identifying novel, disease-associated states within the Microglia and Astrocyte populations and correlating their expression signatures to specific neuroanatomical regions via spatial mapping."*

---

## What We're Using (Data & Tools)

### Data Sources
- **[ssREAD](https://bmblx.bmi.osumc.edu/ssread/)** - Single-cell and spatial RNA-seq database for Alzheimer's Disease (1,053 samples, 277 integrated datasets, 7.3 million cells!)

### Computational Tools
- **Primary Language:** R (v4.x)
- **Core Package:** [Seurat](https://satijalab.org/seurat/) (v5.x) - The workhorse for single-cell RNA-seq analysis



## Learning Resources

### Start Here (Background Knowledge)
🎥 **[StatQuest: A Gentle Introduction to RNA-seq](https://youtu.be/tlf6wYJrwKY?si=DchJbr-J9FhrORXV)** (13 min)
- Watch this first! Josh Starmer explains RNA-seq in the clearest way possible. Don't skip this.

📄 **Key Papers to Read:**
- **Mathys et al. 2019** - [Single-cell transcriptomic analysis of Alzheimer's disease](https://www.nature.com/articles/s41586-019-1195-2) (*Nature*)
- **Lau et al. 2020** - [Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease](https://www.pnas.org/doi/10.1073/pnas.2008762117) (*PNAS*)
- **ssREAD Database Paper** - [A single-cell and spatial RNA-seq database for Alzheimer's disease](https://www.nature.com/articles/s41467-024-49133-z) (*Nature Communications*, 2024)

### Coding Tutorials
📘 **[Seurat v5 Integration Tutorial](https://satijalab.org/seurat/articles/integration_introduction)** (REQUIRED)
- This is the exact workflow we're following. Read it once, don't worry about understanding everything. Then read it again and try the code yourself.

📘 **[Seurat PBMC Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)** (Recommended for beginners)
- If you've never used Seurat before, start here. It's a gentler introduction to the basic workflow.

### Additional YouTube Resources
🎥 **[StatQuest: PCA Main Ideas](https://www.youtube.com/watch?v=HMOI_lkzW08)** (5 min)
🎥 **[StatQuest: UMAP Main Ideas](https://www.youtube.com/watch?v=eN0wFzBA4Sc)** (6 min)
🎥 **[StatQuest: Clustering with DBSCAN](https://www.youtube.com/watch?v=RDZUdRSDOok)** (11 min)


## Installation & Setup

### 1. Install R and RStudio
- Download R: https://cran.r-project.org/
- Download RStudio: https://posit.co/download/rstudio-desktop/

### 2. Install Seurat and Dependencies

```r
# Install Seurat v5
install.packages("Seurat")

# Install other required packages
install.packages(c("dplyr", "ggplot2", "patchwork", "cowplot"))

# For GEO data downloads
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

# For Synapse data downloads
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
```

### 3. Create Synapse Account
- Go to https://www.synapse.org/
- Register for free account (needed for Mathys dataset)




## Project Structure

```
Fa25-Project6-AD-Transcriptomics/
│
├── README.md                    # This file
├── data/                        # Raw and processed data (gitignored)
│   ├── lau2020/
│   └── mathys2019/
├── code/                        # Analysis scripts
│   ├── 01_data_download.R
│   ├── 02_qc_and_integration.R
│   ├── 03_cell_type_annotation.R
│   └── 04_differential_expression.R
├── figures/                     # Generated plots
├── results/                     # DEG lists, tables, etc.
└── final_report/               # Final deliverables
```

---

## Contact & Questions

**Slack Channel:** `#project-6-rna-seq`
**Meeting Time:** Tuesdays 9-10 PM @ Grimes Hall
**Project Lead:** Bhavna

**No question is a dumb question!** Post all your confusion points in the Slack channel so others can learn too.

---

## Acknowledgments

This project uses data from:
- The Religious Orders Study and Memory and Aging Project (ROSMAP)
- The ssREAD database team
- Multiple AD research consortia
