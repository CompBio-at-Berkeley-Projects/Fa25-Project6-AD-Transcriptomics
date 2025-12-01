# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Normalisation,_Selection_and_Scaling
# Potentially good resource

# Ze Zhang, first author of SCINA paper, sadly passed away a couple of
# years ago, so SCINA has not been updated in a while. However, the main
# algorithm should still be relatively stable. 

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
import("ggplot2")

# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD03501.qsave"

# Unfortunately, data redacted for this one. 
# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00203.qsave"
dest <- download_dataset(ssread_data_url, here("ext"))
dataset <- qs::qread(dest)

sig <- list (marker_genes = c("APOE", "APH1B"))

# ChatGPT marker list sample. 
markers <- list(
	Astrocyte    = c("GFAP","AQP4","SLC1A2","GJA1","ALDH1L1","S100B"),
	Excitatory   = c("CAMK2A","SYN3","RBFOX3"),
	Inhibitory   = c("GAD1","GAD2","ERBB4","NXPH1"),
	Microglia    = c("CX3CR1","C1QB","CSF1R","HLA-DRA","CD68","ITGAX"),
	Oligodendrocyte = c("MBP","MOBP","PLP1","MOG","CNP"),
	OPC          = c("PCDH15","MEGF11"),
	Endothelial  = c("CLDN5","FLT1")
)

# More ChatGPT lol
AD_markers <- list(
   Astrocyte = c("AQP4", "GFAP", "SLC1A2", "GJA1", "APOE"),
   Microglia  = c("CSF1R", "CD68", "CX3CR1", "CD74", "APH1B"),
   Oligodendrocyte = c("MBP", "MOG", "PLP1", "CNP"),
   OPC = c("PDGFRA", "VCAN"),
   ExcitatoryNeuron = c("SLC17A7", "CAMK2A", "STMN2"),
   InhibitoryNeuron = c("GAD1", "GAD2"),
   Endothelial = c("CLDN5", "FLT1", "VWF")
)

# So, for reference, layers are known as:
# count: unnormalized
# data: normalized
# scale.data: variance-stabilized
gene_exp_norm <- GetAssayData(object=dataset, assay="RNA", layer="data")
gene_exp_mat <- as.matrix(gene_exp_norm)

# Using SCINA settings from the SCINA github
result <- SCINA(gene_exp_mat, AD_markers,
				max_iter = 100, convergence_n = 10,
				convergence_rate = 0.999, sensitivity_cutoff = 0.9,
				rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')

#View(result$cell_labels)
#View(result$probabilities)

# theoretically, this would write heatmap to pdf
#pdf(here("heatmap.pdf"), width=10, height=8)
plotheat.SCINA(gene_exp_mat, result, AD_markers)
#dev.off()
