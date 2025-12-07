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
	#Pericyte = c("AMBP", "HIGD1B", "PTH1R"),
	#Neuron = c("GLS", "RBFOX3", "CAMK2A"),
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

gen_seurat_dp <- function(result, output)
{
	cat("Generating DimPlot", output, "...\n")
	DimPlot(result$data, reduction="umap")
	ggsave(output)
}

gen_seurat_fp <- function(result, output)
{
	cat("Generating FeaturePlot", output, "...\n")
	FeaturePlot(result$data, reduction="umap", features=our_markers) 
	ggsave(output)
}

gen_scina_dp <- function(result, output)
{
	cat("Generating SCINA DimPlot", output, "...\n")
	DimPlot(result$data, reduction="umap", group.by="scina_labels")
	ggsave(output)
}

gen_scina_bdp <- function(result, output)
{
	cat("Generating better SCINA DimPlot", output, "...\n")
	DimPlot(result$data, reduction="umap",
		group.by="scina_but_better_labels")
	ggsave(output)
}

gen_scina_mp <- function(data, output)
{
	cat("Generating SCINA DimPlot marker plot", output, "...\n")

	dataset <- data$data
	
	cells <- WhichCells(dataset,
		expression=(scina_but_better_labels == "OurMarkers"))

	DimPlot(dataset, cells.highlight=cells, 
		cols.highlight="red", cols="grey80")
	ggsave(output)
}

# processing
scina_process <- function(dataset, gene_exp_mat, output_dir)
{
	# Exclude all cell types with no detectable markers
	markers_final <- lapply(cell_markers, 
		function(m) m[m %in% rownames(dataset)])	
	markers_final <- markers_final[lengths(markers_final) > 0]
		
	# Using SCINA settings from the SCINA github
	# Profiling using system.time instead of Rprof (may change)
	# This SCINA stuff and assignment for the unknowns still need 
	# some tweaking
	cat("Running SCINA...\n")
	result <- SCINA(gene_exp_mat, markers_final, 
		max_iter = 200, convergence_n = 10, convergence_rate = 0.999, 
		sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, 
		log_file="SCINA.log")

	dataset$scina_labels <- result$cell_labels
	names(dataset$scina_labels) <- colnames(dataset)

	better_labels <- dataset$scina_labels

	cluster_majority <- list()

	for (cluster in unique(dataset@meta.data$seurat_clusters)) {
		cells <- which(dataset$seurat_clusters == cluster)
		labels <- dataset$scina_labels[cells]
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

	dataset$scina_but_better_labels <- better_labels
	sink(paste(output_dir, "/celllabels.txt", sep=""))
	print(result$cell_labels)
	sink()

	sink(paste(output_dir, "/probabilities.txt", sep=""))
	print(result$probabilities)
	sink()

	result <- list(
		data=dataset,
		scina_res=result,
		markers=markers_final
	)

	return(result)
}

process_dataset_qsave <- function(qsave_dest, output_dir)
{
	bp_path <- paste(output_dir, "/barplot.png", sep="")
	hp_path <- paste(output_dir, "/heatplot.png", sep="")

	seurat_dp_path <- paste(output_dir, "/seurat_barplot.png", sep="")
	seurat_fp_path <- paste(output_dir, "/seurat_featureplot.png", sep="")

	scina_dp_path  <- paste(output_dir, "/scina_dimplot.png", sep="")
	scina_bdp_path <- paste(output_dir, "/scina_better_dimplot.png", sep="")
	scina_mp_path  <- paste(output_dir, "/scina_markerplot.png", sep="")

	dataset <- qs::qread(qsave_dest)

	# So, for reference, layers are known as:
	# count: unnormalized
	# data: normalized
	# scale.data: variance-stabilized
	gene_exp_norm <- GetAssayData(object=dataset, assay="RNA", layer="data")
	gene_exp_mat <- as.matrix(gene_exp_norm)
	
	result <- scina_process(dataset, gene_exp_mat, output_dir)
	
	# NOTE(wv): too slow to run for debug purposes
	# gen_hp(result, hp_path)
	gen_bp(result, bp_path)

	gen_seurat_dp(result, seurat_dp_path)
	gen_seurat_fp(result, seurat_fp_path)

	gen_scina_dp(result, scina_dp_path)
	gen_scina_bdp(result, scina_bdp_path)
	gen_scina_mp(result, scina_mp_path)
}

process_dataset_from_url <- function(url)
{
	dataset_name <- file_path_sans_ext(basename(url))
	output_dir <- paste("output/", dataset_name, sep="")

	if (!dir.exists("output")) dir.create("output")
	if (!dir.exists(output_dir)) dir.create(output_dir)

	dest <- download_dataset(url, output_dir)
	process_dataset_qsave(dest, output_dir)
}


process_dataset_qsave_seuratv5 <- function(qsave_dest, output_dir)
{
	bp_path <- paste(output_dir, "/barplot.png", sep="")
	hp_path <- paste(output_dir, "/heatplot.png", sep="")

	seurat_dp_path <- paste(output_dir, "/seurat_barplot.png", sep="")
	seurat_fp_path <- paste(output_dir, "/seurat_featureplot.png", sep="")

	scina_dp_path  <- paste(output_dir, "/scina_dimplot.png", sep="")
	scina_bdp_path <- paste(output_dir, "/scina_better_dimplot.png", sep="")
	scina_mp_path  <- paste(output_dir, "/scina_markerplot.png", sep="")

	dataset <- qs::qread(qsave_dest)
	
	layers <- Layers(dataset, assay="RNA")[grepl("^data", 
		Layers(dataset, assay="RNA"))]

	# Exclude all cell types with no detectable markers
	markers_final <- lapply(cell_markers, 
		function(m) m[m %in% rownames(dataset)])	
	markers_final <- markers_final[lengths(markers_final) > 0]
	
	all_labels <- list()
	all_probs <- list()
	for (l in layers) {
		mat <- GetAssayData(dataset, assay="RNA", layer=l)
		#mat <- mat[intersect(markers_final, rownames(mat)), , drop=FALSE]
		mat_scina <- as.matrix(mat)
		
		cat("Running SCINA on", "layer", l, "...\n")
		result <- SCINA(mat_scina, markers_final, 
			max_iter = 200, convergence_n = 10, convergence_rate = 0.999, 
			sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, 
			log_file="SCINA.log")
		cells <- colnames(mat_scina)
		labels <- result$cell_labels
		names(labels) <- cells

		all_labels[[l]] <- labels 
		all_probs[[l]] <- result$probabilities
	}
	
	scina_labels <- c()	
	for (l in names(all_labels)) {
		labels <- all_labels[[l]]
		scina_labels <- c(scina_labels, labels)
	}
	
	#scina_probs <- do.call(rbind, all_probs)
	#rownames(scina_probs) <- names(scina_labels)
	
	#browser()	
	print("sample")	
	# scina_labels <- unlist(all_labels)
	# scina_labels <- unlist(lapply(all_labels, function(x) {x}, use.names=TRUE))

	scina_labels <- scina_labels[colnames(dataset)]
	dataset$scina_labels <- scina_labels

	result <- list(
		data=dataset,
		scina_res=result,
		markers=markers_final
	)


	#result <- scina_process(merged_mat, output_dir)
	
	# NOTE(wv): too slow to run for debug purposes
	# gen_hp(result, hp_path)
	#gen_bp(result, bp_path)

	gen_seurat_dp(result, seurat_dp_path)
	gen_seurat_fp(result, seurat_fp_path)

	gen_scina_dp(result, scina_dp_path)
	#gen_scina_bdp(result, scina_bdp_path)
	#gen_scina_mp(result, scina_mp_path)
}

################################################################################
# NOTE(wv): ~main code~

# Unfortunately, data redacted for this one. 
# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
disease_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00203.qsave"
control_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00202.qsave"

# process_dataset_from_url(disease_url)
# process_dataset_from_url(control_url)

dataset_name <- "GSE157827_AD"
output_dir <- paste("output/", dataset_name, sep="")
if (!dir.exists("output")) dir.create("output")
if (!dir.exists(output_dir)) dir.create(output_dir)

process_dataset_qsave_seuratv5("ext/GSE157827_merged_seurat.qs", output_dir)

