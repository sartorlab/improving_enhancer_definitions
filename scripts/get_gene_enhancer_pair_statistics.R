suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--gene_enhancer_pairs"), action = "store", type = "character", help = "[Required] Path to a file listing the gene-enhancer interactions (one interaction per line; genes on the left, enhancers on the right)"),
  make_option(c("--prefix"), action = "store", type = "character", help = "[Required] Prefix for the names of the files to be printed"),
  make_option(c("--header"), action = "store_true", default = F, help = "[Optional flag] Indicates there is a header at the top of the gene-enhancer interaction list")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, 
                              description = "\nDescription: Accepts a list of gene-enhancer interactions and outputs relevant plots and statistics to files with the prefix given via 'prefix'.")
opts <- parse_args(option_parser)

# for testing
#opts <- list("/Users/peterorchard/Google Drive/Rotations/sartor/interaction_lists/ENCSR752QCX.interactions.fdr.mango.point_to_point.txt", T)
#names(opts) <- c("gene_enhancer_pairs", "header")

#####  VERIFY THE INPUT  #####
if (length(opts$gene_enhancer_pairs) == 0) {
  stop("--gene_enhancer_pairs is a required argument")
}
if (length(opts$prefix) == 0) {
  stop("--prefix is a required argument")
}



#####  MAIN  #####

gene_enhancer_pairs <- ""
if (!opts$header) {
  gene_enhancer_pairs <- read.table(opts$gene_enhancer_pairs, head = F, as.is = T, sep = '\t')
  colnames(gene_enhancer_pairs) <- c("gene", "enhancer")
} else {
  gene_enhancer_pairs <- read.table(opts$gene_enhancer_pairs, head = T, as.is = T, sep = '\t')
}

# ensure we don't double-count anything
gene_enhancer_pairs <- unique(gene_enhancer_pairs)


# generate some very basic statistics that will simply be written to a file
number_unique_enhancers <- length(unique(gene_enhancer_pairs$enhancer))
number_unique_genes <- length(unique(gene_enhancer_pairs$gene))
number_total_interactions <- nrow(gene_enhancer_pairs)


lines <- c(
  paste("There are", number_unique_genes, "genes assigned to at least one enhancer"),
  paste("There are", number_unique_enhancers, "enhancers assigned to at least one gene"),
  paste("There are", number_total_interactions, "total gene-enhancer pairs")
  )

file_conn <- file(paste(opts$prefix, ".stats.txt", sep=''))
writeLines(lines, file_conn)
close(file_conn)


# generate histograms showing the number of interactions per gene/number per enhancer
number_interactions_per_gene <- table(gene_enhancer_pairs$gene)
number_interactions_per_enhancer <- table(gene_enhancer_pairs$enhancer)


png(filename = paste(opts$prefix, ".interaction_histograms.png", sep = ''), width = 5, height = 3, units = "in", res = 300)
par(mfrow=c(1, 2))
hist(number_interactions_per_gene, main = "interactions / gene", xlab = "# interactions", cex = 0.6)
hist(number_interactions_per_enhancer, main = "interactions / enhancer", xlab = "# interactions", cex = 0.6)
dev.off()
par(mfrow=c(1, 1))

# make a plot for distance between gene and enhancer?




