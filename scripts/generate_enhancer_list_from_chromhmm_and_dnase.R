suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--genes"), action = "store", type = "character", help = "[Required] Path to a file listing the genomic intervals that should be assigned to each gene.  Format is one line per locus, tab-separated, chromosome, start, end, gene_id (no header)"),
  make_option(c("--dnase"), action = "store", type = "character", help = "[Required] Path to a file listing the enhancers.  Format is one line per locus, tab-separated, chromosome, start, end, enhancer_id (no header)"),
  make_option(c("--tissue_threshold_dnase"), action = "store", type = "numeric", default = 1, help = "[Optional] An integer.  If the DNase list is used to create the enhancer list, then only DNase sites supported by this number of experiments are included.  There are ~125 experiments in the ENCODE DNase list (Default: 1)."),
  make_option(c("--chrom_hmm_directory"), action = "store", type = "character", help = "[Required] Path to the directory containing the ENCODE chromHMM files. The file names should end in '*HMM.bed.gz'"),
  make_option(c("--extension"), action = "store", type = "numeric", help = "[Required] Number of base pairs to extend enhancers to"),
  make_option(c("--out"), action = "store", type = "character", help = "[Required] File name for the resulting 'pairs' file")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, 
                              description = "\nDescription: Accepts a list of gene loci, enhancers, and an interaction dataset (e.g., ChIA-PET or Hi-C).  Prints a list of gene-enhancer pairs. Generates this list using one of two methods: either by finding the interactions that link one gene to one enhancer (each end of the interaction overlaps with one of the genes/enhancers; 'point_to_point') or by declaring all gene/enhancer pairs encompassed WITHIN an interaction region to be interacting ('encompassing'; useful when e.g. enhancers/genes are expected to be within the loops rather than at the interaction edges, such as would be expected within a TAD).")
opts <- parse_args(option_parser)


#####  VERIFY THE INPUT  #####
if (length(opts$genes) == 0) {
  stop("--genes is a required argument")
}
if (length(opts$chrom_hmm_directory) == 0) {
  stop("--chrom_hmm_directory is a required argument")
}
if (length(opts$extension) == 0) {
  stop("--extension is a required argument")
}
if (length(opts$out) == 0) {
  stop("--out is a required argument")
}

#####  LOAD THE REST OF THE LIBRARIES  #####
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("hash"))
source("../scripts/create_enhancer_list.R", chdir = T)

#####  MAIN  #####

# read in the genes loci
genes <- read.table(opts$genes, head = F, as.is = T, sep = '\t')
colnames(genes) <- c("chromosome", "start", "end", "gene_id")

# load in the chromhmm data for all the tissues
chrom_hmm_files <- list.files(opts$chrom_hmm_directory, pattern = "HMM.bed.gz")
chrom_hmm <- ""
for(file in chrom_hmm_files) {
  cell_type <- gsub("wgEncodeBroadHmm(.*)HMM.bed.gz", "\\1", file, perl = T)
  dat <- read.table(gzfile(paste(opts$chrom_hmm_directory, file, sep='/')), head = F, as.is = T, sep = '\t')
  colnames(dat) <- c("chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
  dat <- dat[,c("chromosome", "start", "end", "name")]
  dat$cell_type <- cell_type
  if (is.data.frame(chrom_hmm)) {
    chrom_hmm <- rbind(chrom_hmm, dat)
  } else {
    chrom_hmm <- dat
  }
}

chrom_hmm <- chrom_hmm[grep("Enhancer", chrom_hmm$name),]
chrom_hmm <- chrom_hmm[,c("chromosome", "start", "end")]

if (length(opts$dnase) > 0) {
  # load in the dnase data
  dnase <- read.table(gzfile(opts$dnase), head = F, as.is = T, sep = '\t')
  colnames(dnase) <- c("chromosome", "start", "end", "name", "score", "floatScore", "sourceCount", "sourceIds")
  dnase$length <- dnase$end - dnase$start
  dnase <- dnase[dnase$sourceCount>=opts$tissue_threshold_dnase,]
  chrom_hmm <- rbind(chrom_hmm, dnase[,c("chromosome", "start", "end")])
}

chrom_hmm <- as.data.frame(chrom_hmm)
chrom_hmm.ranges <- makeGRangesFromDataFrame(chrom_hmm)
chrom_hmm.merged <- reduce(chrom_hmm.ranges)

chrom_hmm <- create_enhancer_list(chrom_hmm.merged, opts$extension)
chrom_hmm <- as.data.frame(chrom_hmm)
chrom_hmm <- chrom_hmm[,c("seqnames", "start", "end")]
colnames(chrom_hmm) <- c("chromosome", "start", "end")

# ditch things that overlap with genes...
genes.ranges <- makeGRangesFromDataFrame(df = genes, seqnames.field = "chromosome", start.field = "start", end.field = "end")
enhancers.ranges <- makeGRangesFromDataFrame(df = chrom_hmm, seqnames.field = "chromosome", start.field = "start", end.field = "end")
overlaps <- findOverlaps(genes.ranges, enhancers.ranges)
kick_out <- unique(subjectHits(overlaps))
chrom_hmm <- chrom_hmm[-1 * kick_out,]
chrom_hmm$id <- paste(chrom_hmm$chromosome, chrom_hmm$start, chrom_hmm$end, sep = ':')

write.table(chrom_hmm[,c("chromosome", "start", "end")], opts$out, quote = F, sep = "\t", row.names = F, col.names = F)

