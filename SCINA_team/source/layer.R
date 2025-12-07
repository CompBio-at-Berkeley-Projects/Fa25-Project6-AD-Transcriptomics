
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


