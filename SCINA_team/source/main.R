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


# Marker list
AD_markers <- list(
   Astrocyte = c("GFAP", "AQP4", "GJA1", "SLC1A2", "FGFR3", "NKAIN4", "AGT", "PLXNB1", "SLC1A3"),
   Microglia  = c("P2RY12", "CSF1R", "C3", "CX3CR1"),
   Oligodendrocyte = c("OLIG2", "MBP", "MOBP", "PLP1", "MYRF", "MAG"),
   OPC = c("VCAN", "SOX8"),
   ExcitatoryNeuron = c("SLC17A6", "SLC17A7", "SATB2"),
   InhibitoryNeuron = c("GAD1", "GAD2"),
   Endothelial = c("CLDN5", "VWF"),
   Pericyte = c("AMBP", "HIGD1B", "PTH1R"),
   Neuron = c("GLS", "RBFOX3", "CAMK2A"),
   OurMarkers = c("APOE", "APP", "ACE", "GRN", "APH1B")
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
				rm_overlap=TRUE, allow_unknown=FALSE, log_file='SCINA.log')

#View(result$cell_labels)
#View(result$probabilities)

# theoretically, this would write heatmap to pdf
#pdf(here("heatmap.pdf"), width=10, height=8)
plotheat.SCINA(gene_exp_mat, result, AD_markers)
#dev.off()

#make SCINA output into a clean dataframe
cell_labels <- result$cell_labels
prob_mat <- result$probabilities

row_index <- match(cell_labels, rownames(prob_mat))
assigned_probs <- rep(NA, length(cell_labels))
valid <- which(!is.na(row_index))

assigned_probs[valid] <- prob_mat[cbind(row_index[valid], valid)]

result_df <- data.frame(
  cell = colnames(prob_mat),
  label = cell_labels,
  probability = assigned_probs
)

#bar plot below
bar_plot <- ggplot(data = result_df, aes(x = cell_labels)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bar_plot
