sctransform <- TRUE

library(Seurat)
library(SeuratData)
library(patchwork)

library(metap)

library(ggplot2) # for plotting
library(cowplot) # for the theme used by Seurat
theme_set(theme_cowplot()) # set the theme 

# already installed on my machine
# InstallData("ifnb")

ifnb <- LoadData("ifnb")

# Split into two layers: one for control, one for stim
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f=ifnb$stim)
ifnb

if (!sctransform) {
	# initial data processing
	# standard log normalization
	ifnb <- NormalizeData(ifnb)
	ifnb <- FindVariableFeatures(ifnb)
	ifnb <- ScaleData(ifnb)
	ifnb <- RunPCA(ifnb)

	# Without integration
	# ifnb <- FindNeighbors(ifnb, dims=1:30, reduction="pca")
	# ifnb <- FindClusters(ifnb, resolution=2, cluster.name="unintegrated_clusters")

	# ifnb <- RunUMAP(ifnb, dims=1:30, reduction="pca", reduction.name="umap.unintegrated")
	# DimPlot(ifnb, reduction="umap.unintegrated", group.by=c("stim", "seurat_clusters"))

	# standard integration
	ifnb <- IntegrateLayers(object=ifnb, method=CCAIntegration, 
		orig.reduction="pca", new.reduction="integrated.cca",
		verbose=FALSE)
	ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

	ifnb <- FindNeighbors(ifnb, reduction="integrated.cca", dims=1:30)
	ifnb <- FindClusters(ifnb, resolution=1)
	ifnb <- RunUMAP(ifnb, dims=1:30, reduction="integrated.cca")

} else {
	# SCTransform normalization
	ifnb <- SCTransform(ifnb)
	ifnb <- RunPCA(ifnb)
	
	ifnb <- IntegrateLayers(object=ifnb, method=CCAIntegration,
		normalization.method="SCT", verbose=FALSE)
	ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

	ifnb <- FindNeighbors(ifnb, reduction="integrated.dr", dims=1:30)
	ifnb <- FindClusters(ifnb, resolution=0.6)
	ifnb <- RunUMAP(ifnb, dims=1:30, reduction="integrated.dr")
}
# Two separate plots
p_dim <- DimPlot(ifnb, reduction="umap", group.by=c("stim", "seurat_annotations"))

# To visualize side-by-side
# p_dim <- DimPlot(ifnb, reduction="umap", split.by="stim")

Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb,
	ident.1="NK", grouping.var="stim",
	verbose=FALSE)
head(nk.markers)

Idents(ifnb) <- factor(Idents(ifnb),
	levels=c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
		"B Activated", "B", "CD8 T", "NK", "T activated",
		"CD4 Naive T", "CD4 Memory T"))
markers_to_plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7",
	"CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2",
	"S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13",
	"IL3RA", "IGJ", "PRSS57")

p_dot <- DotPlot(ifnb, features = markers.to.plot, 
	cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
# Identifying differential expressed genes across conditions
aggregate_ifnb <- AggregateExpression(ifnb, group.by=c("seurat_annotations", "stim"), return.seurat=TRUE)
genes_to_label = c("ISG15", "LY6E", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

# Scatter plots
p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight=genes_to_label)
p2 <- LabelPoints(plot=p1, points=genes_to_label, repel=TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight=genes_to_label)
p4 <- LabelPoints(plot=p3, points=genes_to_label, repel=TRUE)

p_diff <- p2+p4

if (sctransform) {
	ifnb <- PrepSCTFindMarkers(ifnb)
}
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep="_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1="B_STIM", ident.2="B_CTRL", verbose=FALSE)
head(b.interferon.response, n=15)

# Feature plots
p_feature <- FeaturePlot(ifnb, features=c("CD3D", "GNLY", "IFI6"),
	split.by="stim", max.cutoff=3,
	cols=c("grey", "red"), reduction="umap")

# Violin plots
p_vln <- VlnPlot(ifnb, features=c("LYZ", "ISG15", "CXCL10"),
	split.by="stim", group.by="seurat_annotations",
	pt.size=0, combine=FALSE)

# To show in Rstudio plots panel
print(p_dim)
print(p_dot)
print(p_diff)
print(p_feature)
wrap_plots(plots=p_vln, ncol=1)
