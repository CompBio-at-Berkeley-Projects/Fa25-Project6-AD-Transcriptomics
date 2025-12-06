# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Normalisation,_Selection_and_Scaling
# Potentially good resource

# Ze Zhang, first author of SCINA paper, sadly passed away a couple of
# years ago, so SCINA has not been updated in a while. However, the main
# algorithm should still be relatively stable. 

# Assumes the script is run from the working directory
# and not in the source directory. 

source("source/layer.R")

################################################################################
# NOTE(winston): ~packages~

import("BiocManager")

import("SCINA")
import("qs")

import("Seurat")

import("dplyr")

import("ggplot2")
import("cowplot")



################################################################################
# NOTE(winston): ~global variables~
# Marker list
AD_markers <- list(
	Oligodendrocyte = c("OLIG2", "MBP", "MOBP", "PLP1", "MYRF", "MAG"),
	Astrocyte = c("GFAP", "AQP4", "GJA1", "SLC1A2", "FGFR3", "NKAIN4", 
		"AGT", "PLXNB1", "SLC1A3"),
	OPC = c("VCAN", "SOX8"),
	Endothelial = c("CLDN5", "VWF"),
	Microglia  = c("P2RY12", "CSF1R", "C3", "CX3CR1"),
	ExcitatoryNeuron = c("SLC17A6", "SLC17A7", "SATB2"),
	InhibitoryNeuron = c("GAD1", "GAD2"),
	Pericyte = c("AMBP", "HIGD1B", "PTH1R")
	#Neuron = c("GLS", "RBFOX3", "CAMK2A")
#	OurMarkers = c("APOE", "APP", "ACE", "GRN", "APH1B")
)

our_markers <- c("APOE", "APP", "ACE", "GRN", "APH1B")

################################################################################
# NOTE(winston): ~functions~
download_dataset <- function(url, dir)
{
	# Will change to proper download directory later
	# (for testing purposes only)
	download_dest <- paste(dir, basename(url), sep="/")
	if (!file.exists(download_dest)) {
		print(paste("Downloading dataset", basename(url)))
		download.file(ssread_data_url, download_dest, "curl", timeout=600)
	} else {
		cat("Dataset", basename(url), "already downloaded\n")
	}

	return(download_dest)
}

generate_barplot <- function(scina_res, output_file)
{
	#make SCINA output into a clean dataframe
	cell_labels <- scina_res$cell_labels
	prob_mat <- scina_res$probabilities

	row_index <- match(cell_labels, rownames(prob_mat))
	assigned_probs <- rep(NA, length(cell_labels))
	#valid <- which(!is.na(row_index))

	valid <- 1:length(cell_labels)
	assigned_probs[valid] <- prob_mat[cbind(row_index[valid], valid)]

	result_df <- data.frame(
		cell=colnames(prob_mat),
		label=cell_labels,
		probability=assigned_probs
	)

	cat("Generating bar plot", output_file, "\n")
	bar_plot <- ggplot(data = result_df, 
		aes(x = cell_labels)) + geom_bar() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(output_file)
}

generate_heatplot <- function(scina_res, output_file)
{
	#TODO(winston): allow for more filetypes
	png(output_file)
	cat("Generating heatmap", output_file, "\n")
	plotheat.SCINA(gene_exp_mat, result, AD_markers)
	cat("Heatmap generated\n")
	dev.off()
}






################################################################################
# NOTE(winston): ~main code~

# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD03501.qsave"

# Unfortunately, data redacted for this one. 
# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00203.qsave"
dest <- download_dataset(ssread_data_url, "ext")
dataset <- qs::qread(dest)


# So, for reference, layers are known as:
# count: unnormalized
# data: normalized
# scale.data: variance-stabilized
gene_exp_norm <- GetAssayData(object=dataset, assay="RNA", layer="data")
gene_exp_mat <- as.matrix(gene_exp_norm)

# Using SCINA settings from the SCINA github
# Profiling using system.time instead of Rprof (may change)


# This SCINA stuff and assignment for the unknowns still need some tweaking
cat("Running SCINA...\n")
SCINA_prof <- system.time({
	result <- SCINA(gene_exp_mat, AD_markers, 
		max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
		sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, 
		log_file="SCINA.log")
	# result <- SCINA(gene_exp_mat, AD_markers, allow_unknown=FALSE,
		#log_file="SCINA.log")
})

cat("SCINA complete!\n")


dataset$scina_labels <- result$cell_labels
names(dataset$scina_labels) <- colnames(dataset)

cluster_majority <- dataset@meta.data %>%
	group_by(seurat_clusters) %>%
	filter(scina_labels != "unknown") %>%
	summarize(major_label=names(which.max(table(scina_labels))))
better_labels <- dataset$scina_labels

for (cluster in unique(dataset$seurat_clusters)) {
	majority <- cluster_majority$major_label[cluster_majority$seurat_clusters == cluster]
	unknown <- WhichCells(dataset, 
		expression=(scina_labels=="unknown" & seurat_clusters == cluster))
	better_labels[unknown] <- majority
}

dataset$scina_but_better_labels <- better_labels

# heatplot_prof <- system.time(generate_heatplot(result, "output/AD00203_heatplot.png"))
barplot_prof <- system.time(generate_barplot(result, "output/AD00203_barplot.png"))

cat("SCINA time\n")
print(SCINA_prof)
#cat("heatplot time\n")
#print(heatplot_prof)
cat("barplot time\n")
print(barplot_prof)

sink("output/result_celllabels.txt")
result$cell_labels
sink()

sink("output/result_prob.txt")
result$probabilities
sink()


#levels(Idents(dataset)) <- c("Astrocytes", "Endothelial cells", "Excitatory neurons", 
#	"Microglia", "Oligodendrocytes", "Oligodendrocyte precursor cells")

DimPlot(dataset, reduction="umap")
ggsave("output/seurat_dimplot.png")

DimPlot(dataset, reduction="umap", group.by="scina_labels")
ggsave("output/scina_dimplot.png")

DimPlot(dataset, reduction="umap", group.by="scina_but_better_labels")
ggsave("output/better_scina_dimplot.png")


# TODO(winston): fix it so it uses SCINA labelling instead of fricking 
# Seurat labelling
FeaturePlot(dataset, reduction="umap", features=our_markers) 
ggsave("output/seurat_featureplot.png")

