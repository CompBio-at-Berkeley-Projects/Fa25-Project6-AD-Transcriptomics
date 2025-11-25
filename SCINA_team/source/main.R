# Necessary packages
import <- function(package)
{
	mirror_us <- "http://cran.us.r-project.org"
	if(!requireNamespace(package, quietly=TRUE, character.only=TRUE)) {
		install.packages(package, repos=mirror_us)

		# There's gotta be a better way to do this.
		if (package == "BiocManager") {
			BiocManager::install(version = "3.22")
			BiocManager::install("preprocessCore")
		}
	}
	library(package, character.only=TRUE)
}

import("BiocManager")
import("SCINA")
import("qs")

import("Seurat")

download_dataset <- function(url)
{
	# Will change to proper download directory later
	# (for testing purposes only)
	download_dest <- paste("./ext/", basename(url), sep="")
	if (!file.exists(download_dest)) {
		print(paste("Downloading dataset", basename(url)))
		download.file(ssread_data_url, download_dest, "curl", timeout=600)
	} else {
		print(paste("Dataset", basename(url), "already downloaded"))
	}

	return(download_dest)
}

# ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD03501.qsave"
ssread_data_url <- "https://bmblx.bmi.osumc.edu/ssread_download/scrnaseq_qsave/AD00102.qsave"
dest <- download_dataset(ssread_data_url)

dataset <- qs::qread(dest)
