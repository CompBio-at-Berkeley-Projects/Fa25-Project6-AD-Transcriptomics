# Necessary packages
if (!require("BiocManager", quietly=TRUE))
	install.packages("BiocManager", 
					 repos="http://cran.us.r-project.org")
BiocManager::install(version = "3.22")
BiocManager::install("preprocessCore")

if (!require("SCINA", quietly=TRUE))
	install.packages("SCINA", repos="http://cran.us.r-project.org")

library("SCINA")
library("preprocessCore")

# Loading .Rdata different from loading .csv
# Implementing .csv for a more universal format, though I really
# doubt it's standard at all.

gene_exp <- read.csv("ext/example_expmat.csv",
	header=TRUE, row.names=1, stringsAsFactors=FALSE)
sig <- preprocess.signatures("ext/example_signatures.csv")

gene_exp <- log(gene_exp+1)

# Conversion from data frame to matrix
# There must be a more elegant way of doing this, somehow to 
# extract the numerical matrix data out of a data frame. 
gene_exp[] <- lapply(gene_exp, as.numeric)
gene_exp_mat <- as.matrix(gene_exp)

# A lot of the actual data we will be using from ssREAD will be 
# normalized anyways, so hopefully all this will be redundant.
gene_exp_norm <- normalize.quantiles(gene_exp_mat)
# So apparently SCINA needs this to run. Might be the outdated
# parts of the algorithm (last updated 4 years ago) that might be
# necessitating all this work
rownames(gene_exp_norm) <- rownames(gene_exp_mat)
colnames(gene_exp_norm) <- colnames(gene_exp_mat)

# Run default SCINA
res <- SCINA(gene_exp_norm, sig,
	max_iter = 100, convergence_n = 10, 
	convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
	rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')

View(res$cell_labels)
View(res$probabilities)

plotheat.SCINA(gene_exp_norm, res, sig)
