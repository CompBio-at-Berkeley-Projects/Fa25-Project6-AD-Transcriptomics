# http://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Normalisation,_Selection_and_Scaling
# Potentially good resource

# Ze Zhang, first author of SCINA paper, sadly passed away a couple of
# years ago, so SCINA has not been updated in a while. However, the main
# algorithm should still be relatively stable. 

# Assumes the script is run from the working directory
# and not in the source directory. 

source("source/layer.R")

################################################################################
# ~packages~

import("BiocManager")

import("SCINA")
import("qs")

import("Seurat")

import("dplyr")

import("ggplot2")
import("cowplot")

import("tools")
import("Matrix")

################################################################################
# ~global variables~
# Marker list
cell_markers <- list(
	Oligodendrocyte = c("OLIG2", "MBP", "MOBP", "PLP1", "MYRF", "MAG"),
	Astrocyte = c("GFAP", "AQP4", "GJA1", "SLC1A2", "FGFR3", "NKAIN4", 
		"AGT", "PLXNB1", "SLC1A3"),
	OPC = c("VCAN", "SOX8"),
	Endothelial = c("CLDN5", "VWF"),
	Microglia  = c("P2RY12", "CSF1R", "C3", "CX3CR1"),
	ExcitatoryNeuron = c("SLC17A6", "SLC17A7", "SATB2"),
	InhibitoryNeuron = c("GAD1", "GAD2"),
	Pericyte = c("AMBP", "HIGD1B", "PTH1R"),
	Neuron = c("GLS", "RBFOX3", "CAMK2A"),
	OurMarkers = c("APOE", "APP", "ACE", "GRN", "APH1B")
)

our_markers <- cell_markers[["OurMarkers"]] 

################################################################################
# ~functions~
download_dataset <- function(url, dir)
{
	# Will change to proper download directory later
	# (for testing purposes only)
	download_dest <- paste(dir, basename(url), sep="/")
	if (!file.exists(download_dest)) {
		cat("Downloading dataset", basename(url), "\n")
		download.file(url, download_dest, "curl", timeou=600)
	} else {
		cat("Dataset", basename(url), "already downloaded\n")
	}

	return(download_dest)
}

# Building plots
# TODO(wv): Consider whether this should include Unknown or no. 
gen_bp <- function(data, output_file)
{
	#make SCINA output into a clean dataframe
	cell_labels <- data$scina_res$cell_labels
	prob_mat <- data$scina_res$probabilities

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

gen_hp <- function(data, output_file, markers)
{
	#TODO(wv): allow for more filetypes
	png(output_file)
	cat("Generating heatmap...", output_file, "\n")
	plotheat.SCINA(gene_exp_mat, data$scina_res, data$markers)
	dev.off()
}

gen_seurat_dp <- function(dataset, output_dir)
{
	output <- paste(output_dir, "/seurat_dimplot.png", sep="")

	cat("Generating DimPlot", output, "...\n")
	DimPlot(dataset, reduction="umap") + ggtitle("Seurat Clusters")
	ggsave(output)
}

gen_seurat_fp <- function(dataset, output_dir)
{
	all <- paste(output_dir, "/seurat_featureplot_all.png", sep="")

	cat("Generating FeaturePlot", all, "...\n")
	FeaturePlot(dataset, reduction="umap", features=our_markers)
	ggsave(all)

	for (marker in our_markers) {
		output <- paste(output_dir, "/seurat_featureplot_",
			marker, ".png", sep="")
		cat("Generating FeaturePlot", output, "...\n")
		FeaturePlot(dataset, reduction="umap", features=c(marker)) 
		ggsave(output)
	}
}

gen_scina_dp <- function(dataset, output_dir)
{
	output <- paste(output_dir, "/scina_dimplot.png", sep="")
	cat("Generating SCINA DimPlot", output, "...\n")
	all <- DimPlot(dataset, reduction="umap", group.by="scina_labels") +
		ggtitle("Seurat Clusters (SCINA Cell Labelling)")
	ggsave(output)
}

gen_scina_bdp <- function(dataset, output_dir)
{
	output <- paste(output_dir, "/scina_better_dimplot.png", sep="")
	cat("Generating better SCINA DimPlot", output, "...\n")
	DimPlot(dataset, reduction="umap",
		group.by="scina_but_better_labels") +
		ggtitle("Seurat Clusters (SCINA Cell Labelling)")
	ggsave(output)
}

gen_scina_mp <- function(dataset, output_dir)
{
	output  <- paste(output_dir, "/scina_markerplot.png", sep="")
	cat("Generating better SCINA DimPlot marker plot", output, "...\n")

	cells <- WhichCells(dataset,
		expression=(scina_but_better_labels == "OurMarkers"))

	DimPlot(dataset, cells.highlight=cells, 
		cols.highlight="red", cols="grey80") +
		ggtitle("AD Markers Plot")
	ggsave(output)
}

gen_plots <- function(dataset, output_dir)
{
	#bp_path <- paste(output_dir, "/barplot.png", sep="")
	#hp_path <- paste(output_dir, "/heatplot.png", sep="")


	# NOTE(wv): too slow to run for debug purposes
	# gen_hp(result, hp_path)
	# gen_bp(result, bp_path)

	gen_seurat_dp(dataset, output_dir)
	gen_seurat_fp(dataset, output_dir)

	gen_scina_dp(dataset, output_dir)
	gen_scina_bdp(dataset, output_dir)
	gen_scina_mp(dataset, output_dir)
}

# processing
scina_process <- function(dataset, gene_exp_mat)
{
	# Exclude all cell types with no detectable markers
	markers_final <- lapply(cell_markers, 
		function(m) m[m %in% rownames(dataset)])	
	markers_final <- markers_final[lengths(markers_final) > 0]
		
	# Using SCINA settings from the SCINA github
	# Profiling using system.time instead of Rprof (may change)
	# This SCINA stuff and assignment for the unknowns still need 
	# some tweaking
	result <- SCINA(gene_exp_mat, markers_final, 
		max_iter = 200, convergence_n = 10, convergence_rate = 0.999, 
		sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, 
		log_file="SCINA.log")
	
	names(result$cell_labels) <- colnames(gene_exp_mat)

	return(result)
}

scina_process_unknowns <- function(dataset, labels_all)
{
	better_labels <- labels_all
	for (cluster in unique(dataset@meta.data$seurat_clusters)) {
		cells <- which(dataset$seurat_clusters == cluster)
		labels <- labels_all[cells]
		known <- labels[labels != "unknown"]
		
		if (length(known) > 0) {
			tabl <- table(known)
			max <- which.max(tabl)
			majority <- names(max)

			# NOTE(wv): length(tabl) > 1 will consider all unknowns 
			# as "OurMarkers" assuming they are in the same cluster 
			# as cells which are all of "OurMarkers". There might 
			# also be an edge case where there are substantially
			# more of "OurMarkers than there are the other types.
			# Don't know how we should handle that. 

			# The idea right now is that a false negative is
			# preferable, in this case, to a false positive, but we
			# will only see on the long run how practical this will
			# be. 
			if (majority == "OurMarkers" && length(tabl) > 1) {
				new_max <- which.max(tabl[tabl != max(tabl)])
				majority <- names(new_max)
			}
			
			unknown <- cells[labels == "unknown"]
			if (length(unknown) > 0) {
				better_labels[unknown] <- majority
			}
		}
	}

	return(better_labels)
}

process_dataset_qsave <- function(qsave_dest)
{
	cat("Processing dataset at", qsave_dest, "\n")

	dataset_name <- file_path_sans_ext(basename(qsave_dest))
	output_dir <- paste("output/", dataset_name, sep="")

	if (!dir.exists("output")) dir.create("output")
	if (!dir.exists(output_dir)) dir.create(output_dir)

	dataset <- qs::qread(qsave_dest)

	# So, for reference, layers are known as:
	# count: unnormalized
	# data: normalized
	# scale.data: variance-stabilized
	gene_exp_norm <- GetAssayData(object=dataset, assay="RNA", layer="data")
	gene_exp_mat <- as.matrix(gene_exp_norm)
	
	cat("Running SCINA...\n")
	scina_res <- scina_process(dataset, gene_exp_mat)
	dataset$scina_labels <- scina_res$cell_labels
	
	cat("Processing unknowns...\n")
	better_labels <- scina_process_unknowns(dataset, dataset$scina_labels)	
	dataset$scina_but_better_labels <- better_labels

	cat("Generating plots...\n")
	gen_plots(dataset, output_dir)
}

process_dataset_from_url <- function(url, download_dir)
{
	dest <- download_dataset(url, download_dir)
	process_dataset_qsave(dest)
}

process_dataset_qsave_seuratv5 <- function(qsave_dest)
{
	cat("Processing dataset at", qsave_dest, "\n")

	dataset_name <- file_path_sans_ext(basename(qsave_dest))
	output_dir <- paste("output/", dataset_name, sep="")

	if (!dir.exists("output")) dir.create("output")
	if (!dir.exists(output_dir)) dir.create(output_dir)

	dataset <- qs::qread(qsave_dest)
	layers <- Layers(dataset, assay="RNA")[grepl("^data", 
		Layers(dataset, assay="RNA"))]

	# Process dataset by running SCINA on each layer and merging the results when done
	# Technically only works when the cells are different, which they thankfully are
	all_labels <- list()
	all_probs <- list()

	for (l in layers) {
		mat <- GetAssayData(dataset, assay="RNA", layer=l)
		mat_scina <- as.matrix(mat)
		
		cat("Running SCINA on", "layer", l, "...\n")
		scina_res <- scina_process(dataset, mat_scina)

		all_labels[[l]] <- scina_res$cell_labels
		all_probs[[l]] <- scina_res$probabilities
	}
	
	# Merge all labels together	
	scina_labels <- c()	
	for (l in names(all_labels)) {
		labels <- all_labels[[l]]
		scina_labels <- c(scina_labels, labels)
	}

	scina_labels <- scina_labels[colnames(dataset)]
	dataset$scina_labels <- scina_labels
	
	cat("Processing unknowns...\n")
	better_labels <- scina_process_unknowns(dataset, dataset$scina_labels)
	dataset$scina_but_better_labels <- better_labels

	# This theoretically brings the OurMarkers to the top
	# dataset$scina_labels <- factor(dataset$scina_labels,
		#levels=c(setdiff(levels(dataset$scina_labels), "OurMarkers"), "OurMarkers"))
	# dataset$scina_but_better_labels <- factor(dataset$scina_but_better_labels,
		#levels=c(setdiff(levels(dataset$scina_but_better_labels), "OurMarkers"), "OurMarkers"))
	
	gen_plots(dataset, output_dir)
}

################################################################################
# NOTE(wv): ~main code~

# Unfortunately, data redacted for this one. 
# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
disease_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00203.qsave"
control_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00202.qsave"

# process_dataset_from_url(disease_url, "ext")
# process_dataset_from_url(control_url, "ext")

process_dataset_qsave_seuratv5("ext/GSE157827_AD.qs")
process_dataset_qsave_seuratv5("ext/GSE157827_NC.qs")
