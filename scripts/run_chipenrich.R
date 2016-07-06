# take care of options here
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--chipseq"), action = "store", type="character", help = "[Required] file of chip-seq peaks (may be gzipped).  The first three columns must represent chromosome, start, end (header optional)"), # e.g. output from positional_snp_pwm_overlap.pl
  make_option(c("--locus_definitions"), action = "store", type = "character", help = "[Required] file of locus definitions.  Columns are chrom, start, end, gene_id (header required)"), # e.g. output from brooke's pipeline
  make_option(c("--out"), action = "store", help = "[Required] Prefix for results"),
  make_option(c("--n_cores"), action = "store", type = "numeric", default = 1, help = "[Optional] Number of cores to use (default: 1)")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list=option_list, add_help_option = T, prog = "run_chipenrich.R")
opts <- parse_args(option_parser)

#####  VERIFY OPTIONS  #####
if (length(opts$chipseq)==0) {
  stop("Option --chipseq not provided")
}
if (length(opts$locus_definitions)==0) {
  stop("Option --locus_definitions not provided")
}
if (length(opts$out)==0) {
  stop("Option --out not provided")
}


# load up additional necessary libraries
suppressPackageStartupMessages(library("chipenrich"))

# For testing
#opts <- list("/Users/peterorchard/Google Drive/Rotations/sartor/enhancer_definitions/ldef_hg19_5kb_and_enhancers_2kb.bed", "/Users/peterorchard/Google Drive/Rotations/sartor/ENCODE_chipseq/wgEncodeAwgTfbsSydhHuvecMaxUniPk.narrowPeak.gz")
#names(opts) <- c("locus_definitions", "chipseq")
#locusdef_path <- "/Users/peterorchard/Google Drive/Rotations/sartor/enhancer_definitions/ldef_hg19_5kb_and_enhancers_2kb.bed"
#chipseq_peaks <- "/Users/peterorchard/Google Drive/Rotations/sartor/ENCODE_chipseq/wgEncodeAwgTfbsSydhHuvecMaxUniPk.narrowPeak.gz"


### Read in the chipseq data

# check for a header in chipseq file...
chipseq_peaks_first_line <- read.table(gzfile(opts$chipseq), head = F, as.is = T, sep = "\t", nrows = 1)
has_header <- !is.numeric(chipseq_peaks_first_line$V2)

# is the chipseq file gzipped?
is_gzipped <- grep(".gz$", opts$chipseq)

chipseq_peaks <- ""
if (is_gzipped) {
  # no header present
  chipseq_peaks <- read.table(gzfile(opts$chipseq), head = has_header, as.is = T, sep = "\t")
} else {
  chipseq_peaks <- read.table(opts$chipseq, head = has_header, as.is = T, sep = "\t")
}

chipseq_peaks <- chipseq_peaks[,1:3]
colnames(chipseq_peaks) <- c("chrom", "start", "end")

### run chip enrich
results <- chipenrich(chipseq_peaks, out_name = opts$out, locusdef = opts$locus_definitions, n_cores = opts$n_cores)