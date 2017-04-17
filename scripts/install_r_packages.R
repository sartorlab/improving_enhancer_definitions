tryCatch({library(devtools)},
	error = function(e){
		install.packages('devtools', repos = 'http://cran.strudio.com')
	}
)

tryCatch({library(optparse)},
	error = function(e){
		install.packages('optparse', repos = 'http://cran.strudio.com')
	}
)

tryCatch({library(ggplot2)},
	error = function(e){
		install.packages('ggplot2', repos = 'http://cran.strudio.com')
	}
)

tryCatch({library(plyr)},
	error = function(e){
		install.packages('plyr', repos = 'http://cran.strudio.com')
	}
)

tryCatch({library(tidyr)},
	error = function(e){
		install.packages('tidyr', repos = 'http://cran.strudio.com')
	}
)

tryCatch({library(GenomicRanges)},
	error = function(e){
		source("https://bioconductor.org/biocLite.R")
		biocLite("GenomicRanges")
	}
)

tryCatch({library(GenomicRanges)},
	error = function(e){
		source("https://bioconductor.org/biocLite.R")
		biocLite("GenomicRanges")
	}
)

tryCatch({library(BSgenome.Hsapiens.UCSC.hg19)},
	error = function(e){
		source("https://bioconductor.org/biocLite.R")
		biocLite("BSgenome.Hsapiens.UCSC.hg19")
	}
)

tryCatch({library(chipenrich.data)},
	error = function(e){
		devtools::install_github('sartorlab/chipenrich.data')
	}
)

tryCatch({library(chipenrich)},
	error = function(e){
		devtools::install_github('sartorlab/chipenrich')
	}
)
