suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--locus_list"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that are accounted for (tab-separated, 4 columns: chromosome, start, end, class (gene/enhancer))"),
  make_option(c("--peak_list"), action = "store", type = "character", help = "[Required] Path to a file listing the chip-seq peaks (similar format as ENCODE)"),
  make_option(c("--out"), action = "store", type = "character", help = "[Required]")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, 
                              description = "\nDescription: Accepts a list of gene loci, enhancers, and an interaction dataset (e.g., ChIA-PET or Hi-C).  Prints a list of gene-enhancer pairs. Generates this list using one of two methods: either by finding the interactions that link one gene to one enhancer (each end of the interaction overlaps with one of the genes/enhancers; 'point_to_point') or by declaring all gene/enhancer pairs encompassed WITHIN an interaction region to be interacting ('encompassing'; useful when e.g. enhancers/genes are expected to be within the loops rather than at the interaction edges, such as would be expected within a TAD).")
opts <- parse_args(option_parser)


#####  VERIFY THE INPUT  #####
if (length(opts$locus_list) == 0) {
  stop("--locus_list is a required argument")
}
if (length(opts$peak_list) == 0) {
  stop("--peak_list is a required argument")
}
if (length(opts$out) == 0) {
  stop("--out is a required argument")
}


#####  LOAD THE REST OF THE LIBRARIES  #####
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))


#####  MAIN  #####

# read in the loci
#loci <- read.table("/Users/peterorchard/Google Drive/Rotations/sartor/chromhmm/covered_loci.txt", head = F, sep = '\t', as.is = T)
loci <- read.table(opts$locus_list, head = F, as.is = T, sep = '\t')
colnames(loci) <- c("chromosome", "start", "end", "type")
loci <- unique(loci)
total_number_base_pairs <- sum(loci$end - loci$start)


# read in the peak list
peaks <- ""
if (grep(".gz$", "wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak.gz")) {
  peaks <- read.table(gzfile(opts$peak_list), head = F, as.is = T, sep = '\t')
} else {
  peaks <- read.table(opts$peak_list, head = F, as.is = T, sep = '\t')
}

colnames(peaks) <- c("chromosome", "start", "end", "strand", "x1", "x2", "x3", "x4", "x5", "x6")
peaks <- unique(peaks[,c("chromosome", "start", "end")])
# peaks <- read.table(gzfile("/Users/peterorchard/Google Drive/Rotations/sartor/ENCODE_chipseq/wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak.gz"), head = F, as.is = T, sep = '\t')

# now determine the overlaps.  Figure out if each peak overlaps with a gene/enhancer (give preference to gene) or is missed altogether
loci.ranges <- makeGRangesFromDataFrame(df = loci)
peaks.ranges <- makeGRangesFromDataFrame(df = peaks)
overlaps <- findOverlaps(loci.ranges, peaks.ranges)
overlaps <- as.data.frame(cbind(loci[queryHits(overlaps),], peaks[subjectHits(overlaps),]))
overlaps <- overlaps[,c(4, 5, 6, 7)]
overlaps$peak <- paste(overlaps$chromosome, overlaps$start, overlaps$end, sep = ':')

gene_overlaps <- unique(overlaps$peak[overlaps$type=="gene"])
enhancer_overlaps <- unique(overlaps$peak[overlaps$type=="enhancer"])
enhancer_overlaps <- enhancer_overlaps[!enhancer_overlaps %in% gene_overlaps]

number_overlapping_gene <- length(gene_overlaps)
number_overlapping_enhancer <- length(enhancer_overlaps)
number_missed <- nrow(peaks) - number_overlapping_gene - number_overlapping_enhancer
peaks_per_kb <- (number_overlapping_gene + number_overlapping_enhancer) * 1000 / total_number_base_pairs

lines <- c(
    paste("number_overlapping_gene", number_overlapping_gene, sep = '\t'),
    paste("number_overlapping_enhancer", number_overlapping_enhancer, sep = '\t'),
    paste("number_missed", number_missed, sep = '\t'),
    paste("peaks_per_kb", peaks_per_kb, sep = '\t'),
    paste("total_genome_coverage_in_bp", total_number_base_pairs, sep = '\t')
  )

file_conn <- file(opts$out)
writeLines(lines, con = file_conn)
close(file_conn)
