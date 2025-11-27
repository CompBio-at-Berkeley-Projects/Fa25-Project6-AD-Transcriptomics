
# For cool util functions
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

download_dataset <- function(url, dir)
{
	# Will change to proper download directory later
	# (for testing purposes only)
	download_dest <- paste(dir, basename(url), sep="/")
	if (!file.exists(download_dest)) {
		print(paste("Downloading dataset", basename(url)))
		download.file(ssread_data_url, download_dest, "curl", timeout=600)
	} else {
		# print(paste("Dataset", basename(url), "already downloaded"))
	}

	return(download_dest)
}
