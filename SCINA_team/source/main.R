# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Normalisation,_Selection_and_Scaling
# Potentially good resource

# Assumes the script is run from the working directory
# and not in the source directory. 

# This here thing is a little weird to use. 
# TODO(winston): Implement a better method than this
library(here)
source(here("source", "layer.R"))

import("BiocManager")
import("SCINA")
import("qs")

import("Seurat")

# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD03501.qsave"

# Unfortunately, data redacted for this one. 
# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00203.qsave"
dest <- download_dataset(ssread_data_url, here("ext"))
dataset <- qs::qread(dest)

sig <- list (marker_genes = c("APOE", "APH1B"))

# So, for reference, layers are known as:
# count: unnormalized
# data: normalized
# scale.data: variance-stabilized
gene_exp_normalized <- GetAssayData(object=dataset, assay="RNA", layer="data")
gene_exp_mat <- as.matrix(gene_exp_normalized)

result <- SCINA(gene_exp_mat, sig)

#View(result$cell_labels)
#View(result$probabilities)

plotheat.SCINA(gene_exp_normalized, result, sig)

