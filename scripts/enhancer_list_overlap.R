# take care of options here
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--base_enhancer_file"), action = "store", type="character", help = "[Required] Path to file of old enhancers (chrom, start, end; no header)"),
  make_option(c("--new_enhancer_files"), action = "store", type="character", help = "[Required] Comma-separated list of paths to files of new enhancers (chrom, start, end; no header)"),
  make_option(c("--new_enhancer_labels"), action = "store", type="character", help = "[Required] Comma-separated list of labels for each of the new enhancer files (labels for plotting)"),
  make_option(c("--out"), action = "store", help = "[Required] Name of the outfile")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list=option_list, add_help_option = T, prog = "enhancer_list_overlap.R")
opts <- parse_args(option_parser)


#####  VERIFY OPTIONS  #####
if (length(opts$base_enhancer_file)==0) {
  stop("Option --base_enhancer_file not provided")
}
if (length(opts$new_enhancer_files)==0) {
  stop("Option --new_enhancer_files not provided")
}
if (length(opts$new_enhancer_labels)==0) {
  stop("Option --new_enhancer_labels not provided")
}
if (length(opts$out)==0) {
  stop("Option --out not provided")
}

new_enhancer_files <- strsplit(opts$new_enhancer_files, ",")[[1]]
labels <- strsplit(opts$new_enhancer_labels, ",")[[1]]

if (length(new_enhancer_files)!=length(labels)) {
  stop("Number of items for --new_enhancer_files must equal the number of items for --labels")
}


##### FUNCTION DECLARATIONS #####

percentage_overlapping <- function(granges_1, granges_2) {
  ### calculates the percentage of items in granges_1 that show at least 1bp overlap with at least one item in granges_2
  # make sure that the ranges are reduced
  granges_1 <- reduce(granges_1)
  granges_2 <- reduce(granges_2)
  
  # now find the number of ranges in granges_1 that overlap with items in granges_2
  overlaps <- findOverlaps(query = granges_1, subject = granges_2)
  percentage_from_granges_1_overlapping_granges_2 <- length(unique(queryHits(overlaps))) / length(granges_1)
  return(percentage_from_granges_1_overlapping_granges_2)
}

percentage_bp_overlapping <- function(granges_1, granges_2) {
  ### calculates the percentage of bps in granges_1 that are accounted for in granges_2
  # make sure that the ranges are reduced
  granges_1 <- reduce(granges_1)
  granges_2 <- reduce(granges_2)
  
  # first, calculate the total number of bps in granges_1
  granges_1_bps <- sum(end(granges_1) - start(granges_1))
  
  # find the granges_1 bps that are NOT included in granges_2
  in_granges_1_but_not_2 <- setdiff(granges_1, granges_2)
  unshared_granges_1_bps <- sum(end(in_granges_1_but_not_2) - start(in_granges_1_but_not_2))
  
  # final calculation
  percentage_bps_covered <- (granges_1_bps - unshared_granges_1_bps) / granges_1_bps
  
  return(percentage_bps_covered)
}




#####  MAIN  #####

### Plot the overlap between the old enhancer definitions and the new ones
# Measure this in the % of old enhancers with at least one bp overlap with a new enhancer, as well as 
# the % of bps included in the old enhancers that are included in the new enhancers as well
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyr"))


### Load in the old enhancers
old_enhancers <- read.table(opts$base_enhancer_file, head = F, as.is = T, sep = '\t')
colnames(old_enhancers) <- c("chrom", "start", "end")


### Iteratively load in the new enhancers
new_enhancers <- ""
for(i in new_enhancer_files) {
  tmp <- read.table(i, head = F, as.is = T, sep = '\t')
  colnames(tmp) <- c("chrom", "start", "end")
  tmp$label <- labels[match(i, new_enhancer_files)]
  
  if (!is.data.frame(new_enhancers)) {
    new_enhancers <- tmp
  } else {
    new_enhancers <- rbind(new_enhancers, tmp)
  }
}

new_enhancers <- as.data.frame(new_enhancers)

### Now do the overlap calculations

overlaps <- ""
for(i in unique(new_enhancers$label)) {
  granges_1 <- makeGRangesFromDataFrame(df = old_enhancers, seqnames.field = "chrom", start.field = "start", end.field = "end")
  granges_2 <- makeGRangesFromDataFrame(df = new_enhancers[new_enhancers$label==i,], seqnames.field = "chrom", start.field = "start", end.field = "end")
  
  tmp <- data.frame(percent_enhancers_overlapping = c(percentage_overlapping(granges_1, granges_2)), percent_bps_overlapping = c(percentage_bp_overlapping(granges_1, granges_2)))
  tmp$label <- i
  if(is.data.frame(overlaps)) {
    overlaps <- rbind(overlaps, tmp)
  } else {
    overlaps <- tmp
  }
}

# shape the data for ggplot2
overlaps <- tidyr::gather(data = overlaps, unit, percent_overlap, c(1,2))
overlaps$unit[overlaps$unit=="percent_enhancers_overlapping"] <- "enhancers"
overlaps$unit[overlaps$unit=="percent_bps_overlapping"] <- "bps"

# plot
p <- ggplot(overlaps) + geom_bar(aes(x = label, y = percent_overlap, fill = unit), position = "dodge", stat = "identity")
png(filename = opts$out, width = 2 * length(labels), height = 5, units = "in", res = 300)
print(p + ylab("Fraction overlap") + xlab("Enhancer set") + ylim(c(0,1)) + geom_text(aes(x = label, y = percent_overlap, label = round(percent_overlap, 3), group = unit, vjust=-0.2), position=position_dodge(width = 0.9), size = 6))
dev.off()