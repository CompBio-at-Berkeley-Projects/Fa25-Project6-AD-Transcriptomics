# HighRes-AD-Transcriptomics: Unraveling Cellular and Spatial Transcriptomic Heterogeneity in Human Alzheimer's Disease

**Fall 2025 | Computational Biology @ Berkeley | Project 6**

---

## What Are We Actually Doing?

We're using publicly available human brain sequencing data to understand what's happening at the cellular level in Alzheimer's Disease (AD). Specifically, we're trying to figure out:

1. **Which genes are turned on or off in different brain cell types when someone has AD?**
2. **How do specific "vulnerable" cell populations (like Microglia and Astrocytes) change during disease progression?**
3. **Where in the brain are these changes happening?** (using spatial transcriptomics)

Think of it like this: if the brain is a city, and Alzheimer's is like infrastructure breakdown, we want to know exactly which buildings (cell types) are affected, what's going wrong inside them (gene expression changes), and which neighborhoods (brain regions) are hit the hardest.

---

## Team Info

**Duration:** 8 weeks
**Team Size:** 6 members
**Meeting Time:** Tuesdays 9-10 PM @ Grimes Hall
**Project Lead:** Bhavna

---

## The Big Picture Goal

By the end of this project, we'll have:

- ✅ A **Final Scientific Report** summarizing our findings
- ✅ A **Final Presentation** (think conference-style poster/talk)
- ✅ A **Fully-Documented GitHub Repository** with all our code and analysis

**Official Goal Statement:**
*"To utilize human single-cell and spatial RNA sequencing data from the ssREAD/GEO/Synapse portals to characterize the cell-type-specific transcriptional profiles in Alzheimer's Disease, with a focus on identifying novel, disease-associated states within the Microglia and Astrocyte populations and correlating their expression signatures to specific neuroanatomical regions via spatial mapping."*

Translation: We're going to download brain sequencing data, use computational tools to find patterns in how cells behave differently in AD vs. healthy brains, and map those changes to specific locations in the brain.

---

## What We're Using (Data & Tools)

### Data Sources
- **[ssREAD](https://bmblx.bmi.osumc.edu/ssread/)** - Single-cell and spatial RNA-seq database for Alzheimer's Disease (1,053 samples, 277 integrated datasets, 7.3 million cells!)
- **[GEO (Gene Expression Omnibus)](https://www.ncbi.nlm.nih.gov/geo/)** - NCBI's public genomics data repository
- **[Synapse](https://www.synapse.org/)** - Data repository commonly used for neuroscience/AD research (ROSMAP consortium)

### Specific Datasets We're Using

| **Dataset** | **Brain Region** | **Samples** | **Accession** | **Why We Chose It** |
|-------------|------------------|-------------|---------------|---------------------|
| **Lau et al. 2020** | Prefrontal Cortex | 12 AD + 9 Control (~170K nuclei) | [GEO: GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827) | Large cell count, excellent metadata, easy download |
| **Mathys et al. 2019** | Prefrontal Cortex | 24 AD + 24 Control (~80K nuclei) | [Synapse: syn18485175](https://www.synapse.org/#!Synapse:syn18485175) | Gold standard dataset, highly cited, well-annotated |

### Computational Tools
- **Primary Language:** R (v4.x)
- **Core Package:** [Seurat](https://satijalab.org/seurat/) (v5.x) - The workhorse for single-cell RNA-seq analysis
- **Secondary Tools:**
  - Python (Scanpy/Squidpy) for alternative analysis
  - DESeq2/MAST for differential expression
  - Harmony for batch correction
  - GSEA/Enrichr for pathway analysis

---

## How Does This Work? (The Pipeline)

### Week 1-2: **Data Integration & QC**
**What:** Download datasets, clean up low-quality cells, and integrate multiple datasets together
**Why:** Different datasets have "batch effects" (technical noise from different labs/sequencing runs). Integration removes this noise so we can compare apples-to-apples.
**Deliverable:** An integrated cell atlas with preliminary cell type labels

### Week 3-4: **Cell Type Annotation & Clustering**
**What:** Use known marker genes to identify cell types (Microglia, Astrocytes, Neurons, etc.) and cluster them into subpopulations
**Why:** Not all Microglia are the same! Some might be in a "disease-activated" state. We need to find these sub-states.
**Deliverable:** High-resolution cell type map with sub-clusters

### Week 5-6: **Differential Expression Analysis**
**What:** Statistical testing to find genes that are significantly different between AD and Control within each cell type
**Why:** This tells us what's mechanistically going wrong in diseased cells (inflammation? metabolism? cell death pathways?)
**Deliverable:** Lists of differentially expressed genes (DEGs) and pathway enrichment results

### Week 7-8: **Spatial Mapping & Final Report**
**What:** Map our findings onto spatial transcriptomics data to see *where* in the brain these changes occur
**Why:** AD doesn't affect all brain regions equally. Spatial context matters!
**Deliverable:** Final report, presentation, and polished GitHub repo

---

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

---

## Weekly Action Items

### Week 1 (Completed)
- ✅ Read background papers
- ✅ Set up GitHub repo
- ✅ First read-through of Seurat integration tutorial
- ✅ Post questions/confusion points in Slack thread
- ✅ Research answers to your own questions

### Week 2 (Current)
- 🔄 Download Lau and Mathys datasets
- 🔄 Create Seurat objects and perform QC filtering
- 🔄 Implement integration workflow (following tutorial)
- 🔄 Generate preliminary UMAPs and QC plots
- 🔄 **Deliverable:** Integrated cell atlas with batch correction verification

### Week 3-4 (Upcoming)
- Cell type annotation using marker genes
- Subset Microglia and Astrocytes for high-resolution clustering
- Identify disease-associated sub-states

---

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

---

## Important Notes

### "I Don't Have Prior Experience..."
**That's totally fine!** This project is designed for people with little-to-no computational biology background. The weekly structure builds up your skills gradually. Just stay on track with the weekly objectives and watch the YouTube videos.

### "I Don't Understand the Tutorial..."
**Expected!** On your first read, you won't understand most of the jargon. That's why we do:
1. **First read** - Identify what you don't understand
2. **Research** - Look up those terms
3. **Second read** - Now it makes more sense
4. **Implementation** - Run the code yourself

### "Why Are We Doing Integration?"
When you combine datasets from different labs/sequencing runs, they have technical differences called "batch effects." Without integration:
- Cells cluster by *dataset* rather than by *cell type*
- You can't tell if differences are real biology or just technical noise

After integration:
- Cells cluster by *cell type* regardless of which dataset they came from
- Microglia from Dataset 1 and Dataset 2 cluster together
- Now you can fairly compare AD vs. Control

---

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

We're standing on the shoulders of giants who made their data publicly available. Let's make the most of it!

---

**Last Updated:** October 28, 2025
**Project Status:** Week 2 - Data Integration Phase
