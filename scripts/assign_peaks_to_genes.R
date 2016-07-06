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
if (length(opts$peaks)==0) {
  stop("Option --peaks not provided")
}
if (length(opts$loci)==0) {
  stop("Option --loci not provided")
}
if (length(opts$out)==0) {
  stop("Option --out not provided")
}

#####  ADDITIONAL PACKAGES  #####
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ggplot2"))


##### Function declarations #####
has_header <- function(x) {
  # Attempts to tell if the file has a header or not
  # Says that the file has a header if the type of the elements in the first row 
  # don't match the types in the elements of the next 5 rows 
  # (naturally this won't always work, but for common bed files etc it will)
  # Argument: x is the path to a file
  # returns TRUE if it appears that the file has a header; returns FALSE if it appears
  # that the file does not have a header
  
  first_five_rows <- read.table(x, head = T, sep = "\t", as.is = T, nrows = 5)
  first_row <- read.table(x, head = F, sep = "\t", as.is = T, nrows = 1)
  
  for(i in 1:ncol(first_five_rows)) {
    if(class(first_five_rows[,i])!=class(first_row[,i])){
      return(TRUE)
    }
  }
  
  return(FALSE)
}



#####  MAIN  #####

# Load in the peaks
peaks <- ""
if (length(grep(".*\\.gz$", opts$peaks))==1) {
  # gzipped
  peaks <- read.table(gzfile(opts$peaks), header = has_header(opts$peaks), sep = "\t", as.is = T)
} else {
  peaks <- read.table(opts$peaks, header = has_header(opts$peaks), sep = "\t", as.is = T)
}

peaks <- peaks[,1:3]
colnames(peaks) <- c("chromosome", "start", "end")

# Load in the loci
loci <- read.table(opts$loci, header = has_header(opts$loci), sep = "\t", as.is = T)
colnames(loci) <- c("chromosome", "start", "end", "geneid")

# Assign peaks
loci.granges <- makeGRangesFromDataFrame(df = loci, seqnames.field = "chromosome", start.field = "start", end.field = "end")
peaks.granges <- makeGRangesFromDataFrame(df = peaks, seqnames.field = "chromosome", start.field = "start", end.field = "end")

overlaps <- findOverlaps(loci.granges, peaks.granges)
peaks.overlaps <- peaks[subjectHits(overlaps),]
loci.overlaps <- loci[queryHits(overlaps),]

# format for output
peak_ids <- paste(peaks.overlaps$chromosome, peaks.overlaps$start, peaks.overlaps$end, sep = ':')
geneids <- paste(loci.overlaps$geneid)
df <- unique(data.frame(peakid=peak_ids, geneid=geneids))

write.table(x = df, file = opts$out, quote = F, sep = "\t", row.names = F, col.names = opts$no_header)
