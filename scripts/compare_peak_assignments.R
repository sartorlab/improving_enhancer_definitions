# take care of options here
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--peaks"), action = "store", type="character", help = "[Required] List of chip-seq peaks. Tab-separated, header optional.  First three columns presumed to be chrom, start, end"),
  make_option(c("--loci"), action = "store", type="character", help = "[Required] Locus definitions.  chrom, start, end, geneid.  Tab-separated, header optional"),
  make_option(c("--out"), action = "store", help = "[Required] Name of the outfile. Format: peak_id (chrom:start:end), gene_id"),
  make_option(c("--no_header"), action = "store_false", default = T, help = "[Optional] (Flag) Don't include header in the outfile (Default = include header)")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list=option_list, add_help_option = T, prog = "assign_peaks_to_genes.R")
opts <- parse_args(option_parser)


# For testing:
#opts <- list(
#  "/Users/peterorchard/Google Drive/Rotations/sartor/improving_enhancer_definitions/current_definitions/ldef_hg19_5kb_and_enhancers_2kb.bed",
#  "/Users/peterorchard/Google Drive/Rotations/sartor/ENCODE_chipseq/wgEncodeAwgTfbsSydhHuvecGata2UcdUniPk.narrowPeak.gz"
#  )
#names(opts) <- c("loci", "peaks")



#####  VERIFY OPTIONS  #####
if (length(opts$list_1)==0) {
  stop("Option --list_1 not provided")
}
if (length(opts$list_2)==0) {
  stop("Option --list_2 not provided")
}
if (length(opts$out)==0) {
  stop("Option --out not provided")
}

#####  ADDITIONAL PACKAGES  #####
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("plyr"))


#### Temp: load up the lists...
l <- "/Users/peterorchard/Google Drive/Rotations/sartor/improving_enhancer_definitions/test_out.no_header.txt"
list_1 <- read.table(l, header = F, as.is = T, sep = '\t')
list_2 <- read.table(l, header = F, as.is = T, sep = '\t')
colnames(list_1) <- c("peakid", "geneid")
colnames(list_2) <- c("peakid", "geneid")

list_2$geneid <- sample(list_2$geneid)
list_2 <- list_2[sample(1:nrow(list_2), size = 1000),]


### compare the assignments...
all <- unique(rbind(list_1, list_2))

list_1$list_1 <- T
list_2$list_2 <- T

all <- join(all, list_1)
all <- join(all, list_2)
all$list_1[is.na(all$list_1)] <- F
all$list_2[is.na(all$list_2)] <- F

# sanity check
tmp <- apply(all[,c("list_1", "list_2")], 1, sum)
stopifnot(min(tmp)>0)
stopifnot(max(tmp)<=2)

# now classify each row


assignments <- join(list_1, list_2, by = "peakid", type = "full")
assignments$category <- "same"
assignments$category[is.na(assignments$assignment_2)] <- "in_1_not_2"
assignments$category[is.na(assignments$assignment_1)] <- "in_2_not_1"
assignments$category[!is.na(assignments$assignment_1) & !is.na(assignments$assignment_2) & assignments$assignment_2!=assignments$assignment_1] <- "different"

# plot
p <- ggplot(data = assignments) + geom_bar(aes(category))
p
